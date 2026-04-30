// Find Closest (MPI + OpenMP, Distributed Sample Sort)
//   Demonstrate OpenMP + MPI optimizations for BOTH distance
//   computation and global sorting on large datasets (100M+ rows).

//   - Rank 0:
//       * Reads full CSV via load_csv (common.c common.h).
//       * Scatters Aircraft records to all ranks.
//   - All ranks:
//       * Compute phase-specific metric with OpenMP.
//   - Sorting (the interesting part):
//       * use a distributed sample sort across MPI ranks.
//       * Within each rank, we use your **OpenMP sample sort**
//         (same algorithm that beat merge sort in find_closest_omp.c).
//       * After distributed sort, the global list is sorted across ranks:
//         rank 0 has smallest distances, rank P-1 the largest.
//   - Rank 0:
//       * Prints and writes Top-X nearest aircraft from its local
//         (already globally smallest) segment. Top-X is trivial and
//         only for display.
//   mpirun -np 4 ./find_closest_mpi --phase=3 30.4733 -87.1866 25
//           --data-file input/input-data.csv

#define _POSIX_C_SOURCE 200809L
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <omp.h>
#include <mpi.h>

#define SAMPLE_OVERSAMPLE   4   /* for distributed sample sort sampling */

// ---------------------------------------------------------------
//           MATH HELPERS / DISTANCES (driver side)
// ------------------------------------------------------------- 

// distance_from_half is used by common.c for phase 3 display 
double distance_from_half(double a){
    if (a < 0.0) a = 0.0;
    if (a > 1.0) a = 1.0;
    return 2.0 * EARTH_RADIUS_KM * asin(sqrt(a));
}

// Full haversine used in phase 2 and by common.c 
#pragma omp declare simd
double full_haversine(double lat1,double lon1,double lat2,double lon2){
    double latDiff  = deg2rad(lat2 - lat1);
    double lonDiff  = deg2rad(lon_diff(lon2, lon1));
    double rlat1 = deg2rad(lat1), rlat2 = deg2rad(lat2);
    double a = pow(sin(latDiff/2.0),2.0) +
               cos(rlat1)*cos(rlat2)*pow(sin(lonDiff/2.0),2.0);
    return 2.0 * EARTH_RADIUS_KM * asin(sqrt(a));
}

// Half-haversine "a" used in phase 3 
#pragma omp declare simd
double half_haversine(double lat1,double lon1,double lat2,double lon2){
    double latDiff  = deg2rad(lat2 - lat1);
    double lonDiff  = deg2rad(lon_diff(lon2, lon1));
    double rlat1 = deg2rad(lat1), rlat2 = deg2rad(lat2);
    return pow(sin(latDiff/2.0),2.0) +
           cos(rlat1)*cos(rlat2)*pow(sin(lonDiff/2.0),2.0);
}

// Equirectangular approximation used in phase 4 
static inline double equirectangular(double lat1,double lon1,double lat2,double lon2){
    double lonDiff = lon_diff(lon2, lon1);
    double x = deg2rad(lonDiff) * cos(deg2rad((lat1 + lat2)/2.0));
    double y = deg2rad(lat2 - lat1);
    return EARTH_RADIUS_KM * sqrt(x*x + y*y);
}

// -------------------------------------------------------------
//           SORT HELPERS (index-based merge sort)
//  ------------------------------------------------------------

//  Instead of sorting the entire Aircraft structs, we only store
//  the indicies.  This ensures each entry is tracked properly in 
//  sort and not messed up if distances are equal (deterministic)

static inline int cmp_idx_metric(const Aircraft *R, int a, int b){
    if (R[a].metric < R[b].metric) return -1;
    if (R[a].metric > R[b].metric) return  1;
    return (a < b) ? -1 : (a > b);
}


// Merge two adjacent sorted subarrays of src[] into dst[].

static inline void merge_runs(const Aircraft *R,
                              const int *src, int *dst,
                              int left_start,  
                              int mid_start, 
                              int right_end) 
{
    // Read cursors into the left and right runs in src[] 
    int left_pos  = left_start;  // current position in left  run  [left_start, mid_start)  
    int right_pos = mid_start;   // current position in right run  [mid_start,  right_end) 

    int out_pos   = left_start;  // where we write the next smallest element in [left_start, right_end) 

    // While both runs still have elements remaining, choose the smaller
    while (left_pos < mid_start && right_pos < right_end) {
        int left_idx  = src[left_pos];   // index into R for left candidate  
        int right_idx = src[right_pos];  // index into R for right candidate 

        if (R[left_idx].metric < R[right_idx].metric ||
           (R[left_idx].metric == R[right_idx].metric && left_idx <= right_idx)) {

            dst[out_pos++] = left_idx;   // take from left run
            left_pos++;
        } else {
            dst[out_pos++] = right_idx;  // take from right run 
            right_pos++;
        }
    }

    // Copy any remaining elements from the left run (if any)
    while (left_pos < mid_start) {
        dst[out_pos++] = src[left_pos++];
    }

    // Copy any remaining elements from the right run (if any)
    while (right_pos < right_end) {
        dst[out_pos++] = src[right_pos++];
    }
}


// Stable bottom-up mergesort on the index array idx[].
// For small n (<= 32) we fall back to insertion sort.
// For larger n, we do iterative merges with run_width = 1,2,4,8,...
// and optionally parallelize the merge passes with OpenMP.
// I need to config this better..
 
static void merge_sort_idx(const Aircraft *R, int *idx, int n)
{
    if (n <= 1) return;

    int *buf = (int*)malloc((size_t)n * sizeof(int));
    if (!buf) {
        fprintf(stderr, "OOM in merge_sort_idx buffer\n");
        exit(1);
    }

    int *src = idx;
    int *dst = buf;

    // Insertion sort is simple and has very low overhead for tiny n.
    if (n <= 32) {
        for (int i = 1; i < n; ++i) {
            int key = src[i];
            int j = i - 1;
            while (j >= 0 &&
                   (R[src[j]].metric > R[key].metric ||
                   (R[src[j]].metric == R[key].metric && src[j] > key))) {
                src[j+1] = src[j];
                --j;
            }
            src[j+1] = key;
        }
        free(buf);
        return;
    }

    // Bottom-up merge sort
    // On each pass, merge pairs of runs:
    //  [left_start .. mid_start-1] and [mid_start .. right_end-1] into dst[].
    //  After each pass, run_width doubles: 1 -> 2 -> 4 -> 8 -> ...

    for (int width = 1; width < n; width <<= 1) {

        // If we're NOT already inside an OpenMP parallel region and
        // the array is large, parallelize the outer loop over runs.
        if (!omp_in_parallel() && n > 100000) {
            #pragma omp parallel for schedule(static)
            for (int left = 0; left < n; left += (width << 1)) {
                int mid   = left + width;
                int right = left + (width << 1);
                if (mid > n)   mid   = n;
                if (right > n) right = n;

                if (mid >= right ||
                    (R[src[mid-1]].metric < R[src[mid]].metric) ||
                    (R[src[mid-1]].metric == R[src[mid]].metric &&
                     src[mid-1] <= src[mid])) {
                    for (int i = left; i < right; ++i)
                        dst[i] = src[i];
                    continue;
                }
                merge_runs(R, src, dst, left, mid, right);
            }
        } else {
          // If already parallel - or small - then do it serially
            for (int left = 0; left < n; left += (width << 1)) {
                int mid   = left + width;
                int right = left + (width << 1);
                if (mid > n)   mid   = n;
                if (right > n) right = n;

                if (mid >= right ||
                    (R[src[mid-1]].metric < R[src[mid]].metric) ||
                    (R[src[mid-1]].metric == R[src[mid]].metric &&
                     src[mid-1] <= src[mid])) {
                    for (int i = left; i < right; ++i)
                        dst[i] = src[i];
                    continue;
                }
                merge_runs(R, src, dst, left, mid, right);
            }
        }

        //  After this pass, dst has the newly mergd runs
        //  Swap roles so that on the next pass we read from dst and write into src
        int *tmp = src;
        src = dst;
        dst = tmp;
    }

    if (src != idx) {
        memcpy(idx, src, (size_t)n * sizeof(int));
    }
    free(buf);
}

// -----------------------------------------------------------------
//      LOCAL SAMPLE SORT (OpenMP) – used per rank
//  ------------------------------------------------------------- 

// Map a key to a bucket index, given sorted splitters.
// Buckets are [ -inf, split[0] ], (split[0], split[1]], ..., (split[B-2], +inf).
static inline int bucket_of(double key, const double *split, int B){
    int lo = 0, hi = B - 2, where = B - 1;
    while (lo <= hi){
        int mid = (lo + hi) >> 1;
        if (key <= split[mid]) { where = mid; hi = mid - 1; }
        else                   { lo = mid + 1; }
    }
    return where;
}

// OpenMP sample sort over an index array idx[0..n-1] using R[].metric.
// This is essentially the same algorithm as in find_closest_omp.c:
static void sample_sort_idx(const Aircraft *R, int *idx, int n){
    if (n <= 1) return;

    // Choose number of buckets: ~0.5*sqrt(n), minimum 4.
    const int B = (int)fmax(4.0, floor(sqrt((double)n) * 0.5));
    const int s = 4;       // oversampling factor per bucket 
    int S = B * s;
    if (S > n) S = n;

    // 1) Take S evenly spaced samples from idx[], sort them by metric. 
    int *sample = (int*)malloc((size_t)S * sizeof(int));
    if (!sample) {
        fprintf(stderr,"OOM sample\n");
        exit(1);
    }
    double stride = (double)n / (double)S;
    for (int i = 0; i < S; ++i){
        int pos = (int)floor(i * stride);
        if (pos >= n) pos = n - 1;
        sample[i] = idx[pos];
    }

    merge_sort_idx(R, sample, S);

    // 2) Choose B-1 splitters from sorted samples. 
    double *split = (double*)malloc((size_t)(B - 1) * sizeof(double));
    if (!split) {
        fprintf(stderr,"OOM split\n");
        exit(1);
    }
    for (int b = 1; b < B; ++b){
        int pos = (b * S) / B;
        if (pos >= S) pos = S - 1;
        split[b-1] = R[sample[pos]].metric;
    }
    free(sample);

    // 3) Bucket counts and prefix sums. 
    int *bcounts = (int*)calloc((size_t)B, sizeof(int));
    int *boff    = (int*)malloc((size_t)(B + 1) * sizeof(int));
    int *tmp     = (int*)malloc((size_t)n * sizeof(int));
    if (!bcounts || !boff || !tmp){
        fprintf(stderr,"OOM buckets\n");
        exit(1);
    }

    int T = omp_get_max_threads();
    int *lc = (int*)calloc((size_t)T * (size_t)B, sizeof(int));
    if (!lc){
        fprintf(stderr,"OOM local counts\n");
        exit(1);
    }

    // Parallel bucket counting: each thread has its own local counts.
    #pragma omp parallel
    {
        int tid   = omp_get_thread_num();
        int *lc_t = lc + tid * B;

        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i){
            int b = bucket_of(R[idx[i]].metric, split, B);
            lc_t[b]++;
        }
    }

    // Accumulate per-thread counts into global bucket counts. 
    for (int b = 0; b < B; ++b){
        int sum = 0;
        for (int t = 0; t < T; ++t) sum += lc[t * B + b];
        bcounts[b] = sum;
    }

    boff[0] = 0;
    for (int b = 0; b < B; ++b) boff[b+1] = boff[b] + bcounts[b];

    // Per-thread offsets into tmp[]. 
    int *thread_off = (int*)malloc((size_t)T * (size_t)B * sizeof(int));
    if (!thread_off){
        fprintf(stderr,"OOM thread_off\n");
        exit(1);
    }
    for (int b = 0; b < B; ++b){
        int acc = boff[b];
        for (int t = 0; t < T; ++t){
            int c = lc[t * B + b];
            thread_off[t * B + b] = acc;
            acc += c;
        }
    }

    // 4) Parallel scatter into buckets in tmp[].
    #pragma omp parallel
    {
        int tid    = omp_get_thread_num();
        int *off_t = thread_off + tid * B;
        int *cur   = (int*)calloc((size_t)B, sizeof(int));
        if (!cur){
            fprintf(stderr,"OOM cur\n");
            exit(1);
        }

        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i){
            int b = bucket_of(R[idx[i]].metric, split, B);
            int pos = off_t[b] + cur[b]++;
            tmp[pos] = idx[i];
        }
        free(cur);
    }

    // 5) Sort each bucket by metric using merge_sort_idx.
    #pragma omp parallel for schedule(dynamic,1)
    for (int b = 0; b < B; ++b){
        int l = boff[b], r = boff[b+1];
        if (r - l > 1) {
            merge_sort_idx(R, tmp + l, r - l);
        }
    }

    memcpy(idx, tmp, (size_t)n * sizeof(int));

    free(thread_off);
    free(lc);
    free(tmp);
    free(boff);
    free(bcounts);
    free(split);
}

// -----------------------------------------------------------------
//      LOCAL SORT sample_sort_idx
//  ------------------------------------------------------------- 

// sort local array in-place (physically reorder Aircraft[]).
//   - Build idx[0..n-1].
//   - Run sample_sort_idx(local, idx, n) to sort indices by metric.

static void sort_local_array(Aircraft **local_ptr, int *local_n_ptr){
    Aircraft *local = *local_ptr;
    int       n     = *local_n_ptr;
    if (n <= 1 || !local) return;

    int *idx = (int*)malloc((size_t)n * sizeof(int));
    if (!idx) {
        fprintf(stderr, "OOM idx in sort_local_array\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) idx[i] = i;

    // Use sample sort (faster than pure merge sort on large n).
    sample_sort_idx(local, idx, n);

    // Reorder local[] into a new sorted array. 
    Aircraft *sorted = (Aircraft*)malloc((size_t)n * sizeof(Aircraft));
    if (!sorted) {
        fprintf(stderr, "OOM sorted in sort_local_array\n");
        exit(1);
    }
    for (int i = 0; i < n; ++i) {
        sorted[i] = local[idx[i]];
    }

    free(idx);
    free(local);
    *local_ptr   = sorted;
    *local_n_ptr = n;
}

/// -----------------------------------------------------------------
//       DISTRIBUTED SAMPLE SORT (MPI + OpenMP)
//  ------------------------------------------------------------- 

static int cmp_double(const void *a, const void *b){
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

// Distributed sample sort:
//   1) local sort on each rank (OpenMP sample_sort_idx).
//   2)  each rank selects samples of local keys.
//   3) rank 0 gathers all samples, chooses splitters, broadcasts.
//   4) each rank bucketizes its local sorted data and redistributes via MPI_Alltoallv.
//   5) each rank locally sorts its received bucket (again with sample_sort_idx via sort_local_array).
// Now all keys on rank i <= all keys on rank j for i < j:
// the global list is distributed and sorted across ranks.

static void distributed_sample_sort(Aircraft **local_ptr, int *local_n_ptr,
                                    MPI_Comm comm, int rank, int world_size,
                                    MPI_Datatype MPI_AIRCRAFT)
{
    Aircraft *local = *local_ptr;
    int       n     = *local_n_ptr;

    if (world_size == 1) {
        // Single rank: just sort locally and return. 
        sort_local_array(&local, &n);
        *local_ptr   = local;
        *local_n_ptr = n;
        return;
    }

    // 1) local physical sort with OpenMP sample sort ---- 
    sort_local_array(&local, &n);

    // 2) sample selection & splitter choice ---- 
    int S = 0;
    if (n > 0) {
        int base_samples = SAMPLE_OVERSAMPLE * (world_size - 1);
        if (base_samples > n) base_samples = n;
        S = base_samples;
    }

    double *localSamples = NULL;
    if (S > 0) {
        localSamples = (double*)malloc((size_t)S * sizeof(double));
        if (!localSamples) {
            fprintf(stderr, "Rank %d: OOM localSamples\n", rank);
            MPI_Abort(comm, 1);
        }
        // Pick S positions spread across local sorted array.
        for (int k = 0; k < S; ++k) {
            int pos = (int)(((double)(k+1) * n) / (double)(S+1));
            if (pos < 0)      pos = 0;
            if (pos >= n)     pos = n - 1;
            localSamples[k] = local[pos].metric;
        }
    }

    int *allS = NULL;
    if (rank == 0) {
        allS = (int*)malloc((size_t)world_size * sizeof(int));
        if (!allS) {
            fprintf(stderr, "Root: OOM allS\n");
            MPI_Abort(comm, 1);
        }
    }

    MPI_Gather(&S, 1, MPI_INT,
               allS, 1, MPI_INT,
               0, comm);

    double *allSamples = NULL;
    int    *displsS    = NULL;
    int     totalS     = 0;

    if (rank == 0) {
        displsS = (int*)malloc((size_t)world_size * sizeof(int));
        if (!displsS) {
            fprintf(stderr, "Root: OOM displsS\n");
            MPI_Abort(comm, 1);
        }

        displsS[0] = 0;
        for (int r = 1; r < world_size; ++r) {
            displsS[r] = displsS[r-1] + allS[r-1];
        }
        totalS = displsS[world_size-1] + allS[world_size-1];

        if (totalS > 0) {
            allSamples = (double*)malloc((size_t)totalS * sizeof(double));
            if (!allSamples) {
                fprintf(stderr, "Root: OOM allSamples\n");
                MPI_Abort(comm, 1);
            }
        }
    }

    MPI_Gatherv(localSamples, S, MPI_DOUBLE,
                allSamples, allS, displsS, MPI_DOUBLE,
                0, comm);

    if (localSamples) free(localSamples);

    // Root: choose splitters from all samples. 
    double *splitters = NULL;
    if (rank == 0) {
        splitters = (double*)malloc((size_t)(world_size - 1) * sizeof(double));
        if (!splitters) {
            fprintf(stderr, "Root: OOM splitters\n");
            MPI_Abort(comm, 1);
        }

        if (totalS > 0) {
            qsort(allSamples, (size_t)totalS, sizeof(double), cmp_double);
            for (int i = 1; i < world_size; ++i) {
                int idx = (int)(( (long long)i * (long long)totalS ) / world_size);
                if (idx < 0)          idx = 0;
                if (idx >= totalS)    idx = totalS - 1;
                splitters[i-1] = allSamples[idx];
            }
        } else {
            for (int i = 0; i < world_size - 1; ++i) splitters[i] = 0.0;
        }

        if (allSamples) free(allSamples);
        free(allS);
        free(displsS);
    }

    if (rank != 0) {
        splitters = (double*)malloc((size_t)(world_size - 1) * sizeof(double));
        if (!splitters) {
            fprintf(stderr, "Rank %d: OOM splitters\n", rank);
            MPI_Abort(comm, 1);
        }
    }

    MPI_Bcast(splitters, world_size - 1, MPI_DOUBLE, 0, comm);

    // ---- Step 3: bucketization based on splitters & MPI_Alltoallv ----
    int P = world_size;

    int *sendCounts = (int*)calloc((size_t)P, sizeof(int));
    int *sendDispls = (int*)calloc((size_t)P, sizeof(int));
    if (!sendCounts || !sendDispls) {
        fprintf(stderr, "Rank %d: OOM sendCounts/sendDispls\n", rank);
        MPI_Abort(comm, 1);
    }

    // Count buckets per destination rank. 
    if (n > 0) {
        int b = 0;
        for (int i = 0; i < n; ++i) {
            double key = local[i].metric;
            while (b < P - 1 && key > splitters[b]) b++;
            if (b < 0) b = 0;
            if (b >= P) b = P - 1;
            sendCounts[b]++;
        }
    }

    for (int r = 1; r < P; ++r) {
        sendDispls[r] = sendDispls[r-1] + sendCounts[r-1];
    }

    Aircraft *sendBuf = NULL;
    if (n > 0) {
        sendBuf = (Aircraft*)malloc((size_t)n * sizeof(Aircraft));
        if (!sendBuf) {
            fprintf(stderr, "Rank %d: OOM sendBuf\n", rank);
            MPI_Abort(comm, 1);
        }
        int *cursor = (int*)malloc((size_t)P * sizeof(int));
        if (!cursor) {
            fprintf(stderr, "Rank %d: OOM cursor\n", rank);
            MPI_Abort(comm, 1);
        }
        for (int r = 0; r < P; ++r) cursor[r] = sendDispls[r];

        int b = 0;
        for (int i = 0; i < n; ++i) {
            double key = local[i].metric;
            while (b < P - 1 && key > splitters[b]) b++;
            if (b < 0) b = 0;
            if (b >= P) b = P - 1;
            int pos = cursor[b]++;
            sendBuf[pos] = local[i];
        }
        free(cursor);
    }

    free(splitters);
    free(local);  // old local data now lives in sendBuf 

    int *recvCounts = (int*)calloc((size_t)P, sizeof(int));
    int *recvDispls = (int*)calloc((size_t)P, sizeof(int));
    if (!recvCounts || !recvDispls) {
        fprintf(stderr, "Rank %d: OOM recvCounts/recvDispls\n", rank);
        MPI_Abort(comm, 1);
    }

    MPI_Alltoall(sendCounts, 1, MPI_INT,
                 recvCounts, 1, MPI_INT,
                 comm);

    int totalRecv = 0;
    for (int r = 0; r < P; ++r) {
        recvDispls[r] = totalRecv;
        totalRecv    += recvCounts[r];
    }

    Aircraft *new_local = NULL;
    if (totalRecv > 0) {
        new_local = (Aircraft*)malloc((size_t)totalRecv * sizeof(Aircraft));
        if (!new_local) {
            fprintf(stderr, "Rank %d: OOM new_local\n", rank);
            MPI_Abort(comm, 1);
        }
    }

    MPI_Alltoallv(sendBuf, sendCounts, sendDispls, MPI_AIRCRAFT,
                  new_local, recvCounts, recvDispls, MPI_AIRCRAFT,
                  comm);

    free(sendBuf);
    free(sendCounts);
    free(sendDispls);
    free(recvCounts);
    free(recvDispls);

    // ---- Step 4: final local sort of received bucket (OpenMP sample sort) ---- 
    Aircraft *final_local = new_local;
    int       final_n     = totalRecv;
    sort_local_array(&final_local, &final_n);

    *local_ptr   = final_local;
    *local_n_ptr = final_n;
}

// -----------------------------------------------------------------
//           PHASE DRIVER (compute + distributed sort)
//  ------------------------------------------------------------- 

typedef double (*metric_fn_t)(double ref_lat, double ref_lon,
                              double lat, double lon);


// Phase-specific conversion from metric -> km for display 
static double metric_to_km(int phase, double metric){
    if (phase == 3) {
        return distance_from_half(metric);  /* convert 'a' to km */
    }
    return metric;  /* already km for phases 1,2,4 */
}


// Rank 0: print & write top X *distinct* aircraft from its locally sorted chunk.
// After distributed sort, rank 0 holds the globally smallest keys,
// so taking the first X *unique ICAO24* entries is correct.

static void print_and_write_topX_root(int phase,
                                      const char *phase_label,
                                      const char *csv_name,
                                      const Args *args,
                                      Aircraft *local_sorted,
                                      int local_n)
{
    int X = args->topX;
    if (local_n <= 0 || X <= 0) {
        printf("\n--- %s ---\n(no results)\n", phase_label);
        return;
    }

    // We will find up to X unique aircraft (by ICAO24) in metric order.
    int  *unique_pos  = (int*)malloc((size_t)X * sizeof(int));
    char (*seen_icao)[9] = (char (*)[9])calloc((size_t)X, sizeof(*seen_icao));
    if (!unique_pos || !seen_icao) {
        fprintf(stderr, "OOM in print_and_write_topX_root\n");
        free(unique_pos);
        free(seen_icao);
        return;
    }

    int unique_count = 0;

    // Scan sorted array from smallest metric upward, keep first occurrence per ICAO24.
    for (int i = 0; i < local_n && unique_count < X; ++i) {
        const Aircraft *ac = &local_sorted[i];
        const char *icao = ac->icao24;

        // Normalize NULL to empty string.
        if (!icao) icao = "";

        // Check if we've already seen this ICAO24.
        int is_dup = 0;
        for (int j = 0; j < unique_count; ++j) {
            if (strncmp(seen_icao[j], icao, 8) == 0) {
                is_dup = 1;
                break;
            }
        }
        if (is_dup) continue;

        // New aircraft: remember its ICAO and its position in local_sorted.
        strncpy(seen_icao[unique_count], icao, 8);
        seen_icao[unique_count][8] = '\0';
        unique_pos[unique_count] = i;
        unique_count++;
    }

    if (unique_count == 0) {
        printf("\n--- %s ---\n(no distinct aircraft found)\n", phase_label);
        free(unique_pos);
        free(seen_icao);
        return;
    }

    printf("\n--- %s (Top %d distinct aircraft from globally sorted list) ---\n",
           phase_label, unique_count);
    printf("  %-4s %-8s %-8s %9s %10s %8s %10s\n",
           "Rank","ICAO","Callsign","Lat(deg)","Lon(deg)","Alt","Dist(km)");
    printf("  %-4s %-8s %-10s %9s %10s %8s %10s\n",
           "----","--------","--------","---------","----------","--------","----------");

    // Print using unique positions
    for (int j = 0; j < unique_count; ++j) {
        const Aircraft *ac = &local_sorted[unique_pos[j]];
        double dist_km = metric_to_km(phase, ac->metric);
        printf("  %4d %-8s %-8s %9.4f %10.4f %8.0f %10.4f\n",
               j+1,
               dash_if_empty(ac->icao24),
               dash_if_empty(ac->callsign),
               ac->lat, ac->lon, ac->baroaltitude, dist_km);
    }

    // Write the same unique set to CSV
    FILE *fp = fopen(csv_name, "w");
    if (!fp) {
        perror("CSV write");
        free(unique_pos);
        free(seen_icao);
        return;
    }
    fprintf(fp, "rank,time,icao24,callsign,lat,lon,baroaltitude,dist_km\n");
    for (int j = 0; j < unique_count; ++j) {
        const Aircraft *ac = &local_sorted[unique_pos[j]];
        double dist_km = metric_to_km(phase, ac->metric);
        fprintf(fp,"%d,%s,%s,%s,%.6f,%.6f,%.1f,%.3f\n",
                j+1,
                ac->time,
                ac->icao24,
                ac->callsign,
                ac->lat,
                ac->lon,
                ac->baroaltitude,
                dist_km);
    }
    fclose(fp);
    printf("Top-%d distinct aircraft written to %s\n",
           unique_count, csv_name);

    free(unique_pos);
    free(seen_icao);
}


// Main phase driver:
//   - computes metrics on each rank with OpenMP,
//   - runs distributed sample sort with MPI+OpenMP,
//   - rank 0 prints/writes top-X,
//   - logs separate timings for compute vs sort.

static void run_phase_mpi(int phase,
                          const char *phase_label,
                          const char *log_label_base,
                          const char *csv_name,
                          const Args *args,
                          Aircraft **local_ptr, int *local_n_ptr,
                          int global_n,
                          MPI_Comm comm,
                          int rank, int world_size,
                          MPI_Datatype MPI_AIRCRAFT)
{
    const double ref_lat = args->ref_lat;
    const double ref_lon = args->ref_lon;

    metric_fn_t metric = NULL;
    switch (phase) {
        case 1: metric = vincenty_distance_km;   break;
        case 2: metric = full_haversine;   break;
        case 3: metric = half_haversine;     break;
        case 4: metric = equirectangular;   break;
        default: return;
    }

    Aircraft *local = *local_ptr;
    int       local_n = *local_n_ptr;

    // ---------- 1) COMPUTE METRICS (timed, OpenMP inside each rank) ---------- 

    MPI_Barrier(comm);
    double t0 = MPI_Wtime();

    if (local_n > 0 && metric) {
        #pragma omp parallel for simd schedule(static)
        for (int i = 0; i < local_n; ++i) {
            local[i].metric = metric(ref_lat, ref_lon,
                                     local[i].lat, local[i].lon);
        }
    }

    MPI_Barrier(comm);
    double t1 = MPI_Wtime();

    // ---------- 2) DISTRIBUTED SAMPLE SORT (timed, MPI + OpenMP) ---------- 

    MPI_Barrier(comm);
    double t2 = MPI_Wtime();

    distributed_sample_sort(&local, &local_n, comm, rank, world_size, MPI_AIRCRAFT);

    MPI_Barrier(comm);
    double t3 = MPI_Wtime();

    *local_ptr   = local;
    *local_n_ptr = local_n;

    // ---------- 3) ROOT: print / write Top-X (trivial) ---------- 

    if (rank == 0) {
        print_and_write_topX_root(phase, phase_label, csv_name, args,
                                  local, local_n);
    }

    // ---------- 4) REDUCE TIMINGS & LOG ---------- 

    double local_compute = t1 - t0;          // metrics only 
    double local_sort    = t3 - t2;          // distributed sample sort 
    double max_compute   = 0.0;
    double max_sort      = 0.0;

    MPI_Reduce(&local_compute, &max_compute, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&local_sort,    &max_sort,    1, MPI_DOUBLE, MPI_MAX, 0, comm);

    if (rank == 0) {
        double ms_comp = max_compute * 1e3;
        double ms_sort = max_sort    * 1e3;
        int    threads = omp_get_max_threads();
        int    cores   = world_size * threads;

        printf("[Time] %s COMPUTE (metrics only): %.3f ms across %d ranks x %d threads\n",
               log_label_base, ms_comp, world_size, threads);
        printf("[Time] %s SORT+MPI (distributed sample sort): %.3f ms across %d ranks x %d threads\n",
               log_label_base, ms_sort, world_size, threads);

        char label_comp[128];
        char label_sort[128];
        snprintf(label_comp, sizeof(label_comp), "%s_Compute", log_label_base);
        snprintf(label_sort, sizeof(label_sort), "%s_Sort",    log_label_base);

        log_timing(LOG_FILE, label_comp, global_n, cores, ms_comp);
        log_timing(LOG_FILE, label_sort, global_n, cores, ms_sort);
    }
}

// -----------------------------------------------------------------
//             UTIL: TOPOLOGY & TIMER
//  ------------------------------------------------------------- 

static void print_timer_res_root(int rank){
    if (rank == 0) {
        double tick = omp_get_wtick();
        printf("OpenMP timer resolution: %.3f ns\n", tick * 1e9);
    }
}

static void print_rank_thread_info(MPI_Comm comm, int rank, int world_size){
    MPI_Barrier(comm);
    if (rank == 0) {
        printf("\nMPI/OpenMP topology:\n");
    }
    for (int r = 0; r < world_size; ++r) {
        MPI_Barrier(comm);
        if (rank == r) {
            char name[MPI_MAX_PROCESSOR_NAME];
            int len = 0;
            MPI_Get_processor_name(name, &len);
            int threads = omp_get_max_threads();
            printf("  Rank %d on %s: %d OpenMP threads\n", r, name, threads);
            fflush(stdout);
        }
    }
    MPI_Barrier(comm);
    if (rank == 0) printf("\n");
}

// -------------------------------------------------------------
//             MAIN
// ------------------------------------------------------------- 

int main(int argc, char *argv[])
{
    int rank, world_size;
    int provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (rank == 0) {
        printf("MPI ranks: %d, MPI threading level provided=%d\n",
               world_size, provided);
    }

    Args args = parse_args(argc, argv);
    if (args.phase < 1 || args.phase > 4) {
        if (rank == 0) {
            fprintf(stderr,
                "Usage: %s --phase {1|2|3|4} [lat lon X] [--data-file <path>] [--use-rows N]\n", argv[0]);
            fprintf(stderr,
                "Example: %s --phase=3 30.4733 -87.1866 25 --data-file input/input-data.csv --use-rows 10000000\n",
                argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    // Threads per rank for filename & logging 
    int threads = omp_get_max_threads();

    // Short phase tag for filename 
    const char *phase_tag =
        (args.phase == 1) ? "phase1" :
        (args.phase == 2) ? "phase2" :
        (args.phase == 3) ? "phase3" : "phase4";

    //  Example: output/phase3-output-mpi-4n-8t.csv
    //   where 4n = 4 MPI ranks ("nodes"), 8t = 8 threads per rank 
    static char csv_name[256];
    snprintf(csv_name, sizeof(csv_name),
             "output/%s-output-mpi-%dn-%dt.csv",
             phase_tag, world_size, threads);



    if (rank == 0) {
        if (ensure_dir("output") != 0) {
            perror("ensure_dir(\"output\")");
        }

        printf("Phase: %d\n", args.phase);
        printf("MPI ranks (nodes): %d, threads per rank: %d (total logical cores ~ %d)\n",
               world_size, threads, world_size * threads);
        printf("Ref: (%.4f, %.4f), TopX=%d\n",
               args.ref_lat, args.ref_lon, args.topX);
        printf("Data file: %s\n", args.data_file);
        printf("Output CSV: %s\n", csv_name);
        if (args.use_rows > 0) {
            printf("Use rows (max, global): %d\n", args.use_rows);
        }

    }


    // Rank 0 loads full CSV; others start empty.
    Aircraft *global_all   = NULL;
    int       global_count = 0;
    int       capacity     = 0;
    int       rc           = 0;

    if (rank == 0) {

        // If data_file looks like a glob pattern, use load_csv_pattern.
        if (strpbrk(args.data_file, "*?[]") != NULL) {
            rc = load_csv_pattern(args.data_file,
                                  &global_all, &global_count, &capacity);
        } else {
            rc = load_csv(args.data_file,
                          &global_all, &global_count, &capacity);
        }

        if (rc != 0) {
            fprintf(stderr, "Failed to load from '%s'\n", args.data_file);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        printf("Loaded %d records\n", global_count);

                /* Apply global cap from --use-rows if requested */
        if (args.use_rows > 0 && global_count > args.use_rows) {
            printf("Capping to first %d rows due to --use-rows\n",
                   args.use_rows);
            global_count = args.use_rows;
        }

    }

    // Broadcast global_count so all ranks know work size. 
    MPI_Bcast(&global_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (global_count == 0) {
        if (rank == 0 && global_all) free(global_all);
        MPI_Finalize();
        return 0;
    }

    // Block distribution of rows (almost-even chunk per rank).
    int base = global_count / world_size;
    int rem  = global_count % world_size;
    int local_n = base + ((rank < rem) ? 1 : 0);

    // Scatterv metadata on root. 
    int *sendCounts = NULL;
    int *sendDispls = NULL;
    if (rank == 0) {
        sendCounts = (int*)malloc((size_t)world_size * sizeof(int));
        sendDispls = (int*)malloc((size_t)world_size * sizeof(int));
        if (!sendCounts || !sendDispls) {
            fprintf(stderr, "Root: OOM for Scatterv metadata\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        int offset = 0;
        for (int r = 0; r < world_size; ++r) {
            int ln = base + ((r < rem) ? 1 : 0);
            sendCounts[r] = ln;
            sendDispls[r] = offset;
            offset += ln;
        }
    }

    // MPI datatype for Aircraft (byte-wise contiguous). 
    MPI_Datatype MPI_AIRCRAFT;
    MPI_Type_contiguous((int)sizeof(Aircraft), MPI_BYTE, &MPI_AIRCRAFT);
    MPI_Type_commit(&MPI_AIRCRAFT);

    // Allocate local array and scatter.
    Aircraft *local = NULL;
    if (local_n > 0) {
        local = (Aircraft*)malloc((size_t)local_n * sizeof(Aircraft));
        if (!local) {
            fprintf(stderr, "Rank %d: OOM for local array (n=%d)\n", rank, local_n);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    MPI_Scatterv(global_all, sendCounts, sendDispls, MPI_AIRCRAFT,
                 local, local_n, MPI_AIRCRAFT,
                 0, MPI_COMM_WORLD);

    if (rank == 0) {
        free(global_all);
        free(sendCounts);
        free(sendDispls);
    }

    print_timer_res_root(rank);
    print_rank_thread_info(MPI_COMM_WORLD, rank, world_size);


    // Decide labels for the chosen phase.  csv_name was built earlier.
    const char *phase_label = NULL;
    const char *log_label   = NULL;

    switch (args.phase) {
        case 1:
            phase_label = "Phase 1 (Vincenty, MPI+OMP)";
            log_label   = "Phase1_Vincenty_MPI_OMP";
            break;
        case 2:
            phase_label = "Phase 2 (Full Haversine, MPI+OMP)";
            log_label   = "Phase2_FullHaversine_MPI_OMP";
            break;
        case 3:
            phase_label = "Phase 3 (Half Haversine 'a', MPI+OMP)";
            log_label   = "Phase3_HalfHaversine_MPI_OMP";
            break;
        case 4:
            phase_label = "Phase 4 (Equirectangular, MPI+OMP)";
            log_label   = "Phase4_Equirectangular_MPI_OMP";
            break;
        default:
            break;
    }



    // Run phase: metric compute + distributed sample sort + top-X display.
    run_phase_mpi(args.phase,
                  phase_label,
                  log_label,
                  csv_name,
                  &args,
                  &local, &local_n,
                  global_count,
                  MPI_COMM_WORLD,
                  rank, world_size,
                  MPI_AIRCRAFT);

    if (local) free(local);
    MPI_Type_free(&MPI_AIRCRAFT);
    MPI_Finalize();
    return 0;
}

