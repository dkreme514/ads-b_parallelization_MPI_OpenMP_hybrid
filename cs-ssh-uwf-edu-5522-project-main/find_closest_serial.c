// Find Closest (Serial)
//   - No OpenMP anywhere.
//   - Phases:
//       1: Vincenty distance (km)
//       2: Full haversine (km)
//       3: Half-haversine "a" (0..1)
//       4: Equirectangular distance (km) with adaptive pruning

#define _POSIX_C_SOURCE 200809L
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "common.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>

//   ---------- app config ---------- 

#define EPSILON     0.05
#define ODD_EVEN_MAX_N 150000   // max n to run even–odd sort 

// ---------- math (haversine & helpers) ---------- 


// Convenience wrapper 

// full haversine / half haversine / distance_from_half 
double full_haversine(double lat1,double lon1,double lat2,double lon2){
    double dlat  = deg2rad(lat2 - lat1);
    double dlon  = deg2rad(lon_diff(lon2, lon1));
    double rlat1 = deg2rad(lat1), rlat2 = deg2rad(lat2);
    double a = pow(sin(dlat/2.0),2.0) +
               cos(rlat1)*cos(rlat2)*pow(sin(dlon/2.0),2.0);
    return 2.0 * EARTH_RADIUS_KM * asin(sqrt(a));
}

double half_haversine(double lat1,double lon1,double lat2,double lon2){
    double dlat  = deg2rad(lat2 - lat1);
    double dlon  = deg2rad(lon_diff(lon2, lon1));
    double rlat1 = deg2rad(lat1), rlat2 = deg2rad(lat2);
    return pow(sin(dlat/2.0),2.0) +
           cos(rlat1)*cos(rlat2)*pow(sin(dlon/2.0),2.0);
}

double distance_from_half(double a){
    if (a < 0) a = 0;
    if (a > 1) a = 1;
    return 2.0 * EARTH_RADIUS_KM * asin(sqrt(a));
}

static inline double equirectangular(double lat1,double lon1,double lat2,double lon2){
    double dlon = lon_diff(lon2, lon1);
    double x = deg2rad(dlon) * cos(deg2rad((lat1 + lat2)/2.0));
    double y = deg2rad(lat2 - lat1);
    return EARTH_RADIUS_KM * sqrt(x*x + y*y);
}

// ---------- sorts over index by R[].metric (serial) ----------
static inline int cmp_idx_metric(const Aircraft *R, int a, int b){
    if (R[a].metric < R[b].metric) return -1;
    if (R[a].metric > R[b].metric) return  1;
    return (a<b)?-1:(a>b);
}

static void odd_even_sort_idx(const Aircraft *R, int *idx, int n){
    int sorted = 0;
    while (!sorted){
        sorted = 1;
        for (int i=1; i+1<n; i+=2){
            if (cmp_idx_metric(R, idx[i], idx[i+1]) > 0){
                int t=idx[i]; idx[i]=idx[i+1]; idx[i+1]=t; sorted = 0;
            }
        }
        for (int i=0; i+1<n; i+=2){
            if (cmp_idx_metric(R, idx[i], idx[i+1]) > 0){
                int t=idx[i]; idx[i]=idx[i+1]; idx[i+1]=t; sorted = 0;
            }
        }
    }
}

static void merge_idx(const Aircraft *R, int *idx, int *buf, int l, int m, int r){
    int i=l, j=m, k=l;
    while (i<m && j<r){
        if (cmp_idx_metric(R, idx[i], idx[j]) <= 0) buf[k++] = idx[i++];
        else                                         buf[k++] = idx[j++];
    }
    while (i<m) buf[k++] = idx[i++];
    while (j<r) buf[k++] = idx[j++];
    for (int x=l; x<r; ++x) idx[x] = buf[x];
}
static void mergesort_idx_rec(const Aircraft *R, int *idx, int *buf, int l, int r){
    if (r - l <= 32){
        for (int i=l+1; i<r; ++i){
            int key = idx[i], j = i-1;
            while (j>=l && cmp_idx_metric(R, idx[j], key) > 0){
                idx[j+1]=idx[j]; --j;
            }
            idx[j+1]=key;
        }
        return;
    }
    int m = l + ((r-l)>>1);
    mergesort_idx_rec(R, idx, buf, l, m);
    mergesort_idx_rec(R, idx, buf, m, r);
    if (cmp_idx_metric(R, idx[m-1], idx[m]) <= 0) return;
    merge_idx(R, idx, buf, l, m, r);
}
static void merge_sort_idx(const Aircraft *R, int *idx, int n){
    int *buf = (int*)malloc(sizeof(int)*n);
    if (!buf){ fprintf(stderr,"OOM merge buf\n"); exit(1); }
    mergesort_idx_rec(R, idx, buf, 0, n);
    free(buf);
}

// ---------- reporting ---------- 
static bool same_topX2(const int *A, const int *B, int X){
    for (int i=0;i<X;i++) if (A[i]!=B[i]) return false;
    return true;
}

static void run_two_sorts_fullarray_and_report(Aircraft *R, int count, int X,
                                               int phase, double ref_lat, double ref_lon)
{

    if (count <= 0 || X <= 0){
        puts("\n(no results)");
        return;
    }
    int limit = (X < count) ? X : count;

    // base index [0..count-1] 
    int *base = (int*)malloc(sizeof(int)*count);
    if (!base){
        fprintf(stderr,"OOM idx base\n");
        exit(1);
    }
    for (int i=0; i<count; i++) base[i] = i;

    int *idx_odd   = NULL;
    int *idx_merge = (int*)malloc(sizeof(int)*count);
    if (!idx_merge){
        fprintf(stderr,"OOM idx merge\n");
        free(base);
        exit(1);
    }
    memcpy(idx_merge, base, sizeof(int)*count);

    double t0, t1;

    // Optionally run even–odd sort only if the array is not too large 
    if (count <= ODD_EVEN_MAX_N) {
        idx_odd = (int*)malloc(sizeof(int)*count);
        if (!idx_odd){
            fprintf(stderr,"OOM idx odd\n");
            free(idx_merge);
            free(base);
            exit(1);
        }
        memcpy(idx_odd, base, sizeof(int)*count);

        t0 = wall_time();
        odd_even_sort_idx(R, idx_odd, count);
        t1 = wall_time();
        double ms_odd = (t1 - t0) * 1000.0;
        printf("[Time] Even-Odd Transposition Sort (full array): %.3f ms\n", ms_odd);
    } else {
        printf("[Skip] Even-Odd Transposition Sort for n=%d (> %d)\n",
               count, ODD_EVEN_MAX_N);
    }

    // Always do merge sort 
    t0 = wall_time();
    merge_sort_idx(R, idx_merge, count);
    t1 = wall_time();
    double ms_merge = (t1 - t0) * 1000.0;
    printf("[Time] Merge Sort (full array)                : %.3f ms\n", ms_merge);

    // Compare tops if we have an odd–even result, otherwise just show merge 
    if (idx_odd && same_topX2(idx_odd, idx_merge, limit)){
        print_top_from_idx(R, idx_merge, count, X, phase, ref_lat, ref_lon,
                           "Top-X (both sorts agree)");
    } else {
        if (idx_odd) {
            print_top_from_idx(R, idx_odd,   count, X, phase, ref_lat, ref_lon,
                               "Top-X [Even-Odd]");
        }
        print_top_from_idx(R, idx_merge, count, X, phase, ref_lat, ref_lon,
                           "Top-X [Merge]");
    }

    const char *auto_name =
        (phase==1) ? "output/phase1-output-serial.csv" :
        (phase==2) ? "output/phase2-output-serial.csv" :
        (phase==3) ? "output/phase3-output-serial.csv" :
                     "output/phase4-output-serial.csv";

    // CSV uses idx_merge; write_csv_from_idx will re-apply the unique-ICAO logic 
    write_csv_from_idx(auto_name, R, idx_merge, count, X, phase, ref_lat, ref_lon);

    free(idx_merge);
    free(idx_odd);
    free(base);

}

// ---------- phases ---------- 
// Phase 1: vincenty_distance_km, metric = km 
static void phase1(Aircraft *aircrafts,int count,int X,double ref_lat,double ref_lon){
    double t0 = wall_time();
    for(int i=0;i<count;i++){
        aircrafts[i].metric = vincenty_distance_km(ref_lat, ref_lon,
                                                   aircrafts[i].lat, aircrafts[i].lon);
    }
    double t1 = wall_time();
    double elapsed_ms = (t1 - t0) * 1000.0;
    printf("Phase 1 (Vincenty) metric runtime: %.3f ms\n", elapsed_ms);
    log_timing(LOG_FILE,"Phase1_Vincenty",count,1,elapsed_ms);
    run_two_sorts_fullarray_and_report(aircrafts, count, X, 1, ref_lat, ref_lon);
}

// Phase 2: full haversine, metric = km 
static void phase2(Aircraft *aircrafts,int count,int X,double ref_lat,double ref_lon){
    double t0=wall_time();
    for(int i=0;i<count;i++){
        aircrafts[i].metric = full_haversine(ref_lat,ref_lon,
                                             aircrafts[i].lat,aircrafts[i].lon);
    }
    double t1=wall_time();
    double elapsed_ms = (t1 - t0) * 1000.0;
    printf("Phase 2 (Full Haversine) metric runtime: %.3f ms\n",elapsed_ms);
    log_timing(LOG_FILE,"Phase2_FullHaversine",count,1,elapsed_ms);
    run_two_sorts_fullarray_and_report(aircrafts, count, X, 2, ref_lat, ref_lon);
}

// Phase 3: half_haversine, metric = 'a' (0..1) 
static void phase3(Aircraft *aircrafts,int count,int X,double ref_lat,double ref_lon){
    double t0 = wall_time();
    for(int i=0;i<count;i++){
        aircrafts[i].metric = half_haversine(ref_lat,ref_lon,
                                             aircrafts[i].lat,aircrafts[i].lon);
    }
    double t1 = wall_time();
    double elapsed_ms = (t1 - t0) * 1000.0;
    printf("Phase 3 (Half Haversine 'a') metric runtime: %.3f ms\n", elapsed_ms);
    log_timing(LOG_FILE,"Phase3_HalfHaversine",count,1,elapsed_ms);
    run_two_sorts_fullarray_and_report(aircrafts, count, X, 3, ref_lat, ref_lon);
}

static int cmp_hit_dist(const void *p1,const void *p2){
    double a=((const Hit*)p1)->dist_km, b=((const Hit*)p2)->dist_km;
    return (a>b) - (a<b);
}

// Phase 4: adaptive pruning + equirectangular, metric = equirectangular km 
static void phase4(Aircraft *arr,int count,int X,double ref_lat,double ref_lon){
    double coslat = cos(deg2rad(ref_lat));
    if (fabs(coslat) < 1e-6) coslat = 1e-6;
    int use_lon = fabs(ref_lat) < 80.0;

    if (X <= 0 || count <= 0){ 
	puts("No candidates after streaming."); 
	return; 
    }

    double t0 = wall_time();

    Hit *heap = (Hit*)malloc((size_t)X * sizeof(Hit));
    if (!heap){ fprintf(stderr,"OOM heap\n"); return; }
    int k = 0;
    double cutoff_km = INFINITY;
    double lat_thresh_deg = INFINITY;
    double lon_thresh_deg = use_lon ? INFINITY : 360.0;
    int bbox_survivors = 0;

    for (int i=0; i<count; ++i){
        if (k == X){
            double dlat = fabs(arr[i].lat - ref_lat);
            if (dlat > lat_thresh_deg) continue;
            if (use_lon){
                double dlon = fabs(lon_diff(arr[i].lon, ref_lon));
                if (dlon > lon_thresh_deg) continue;
            }
        }

        double d_eq = equirectangular(ref_lat, ref_lon, arr[i].lat, arr[i].lon);
        ++bbox_survivors;


        // --- check if this ICAO is already in the heap --- 
        int dup_idx = -1;
        for (int h = 0; h < k; ++h){
            if (strcmp(arr[heap[h].idx].icao24, arr[i].icao24) == 0){
                dup_idx = h;
                break;
            }
        }

        if (dup_idx >= 0){
            // already have this aircraft; keep only the closer one 
            if (d_eq < heap[dup_idx].dist_km){
                heap[dup_idx].dist_km = d_eq;
                // distance decreased -> fix heap by sifting down 
                heap_sift_down(heap, k, dup_idx);
                if (k == X){
                    cutoff_km      = heap_top(heap).dist_km * (1.0 + EPSILON);
                    lat_thresh_deg = cutoff_km / 111.0;
                    lon_thresh_deg = use_lon ? cutoff_km / (111.0 * coslat) : 360.0;
                }
            }
            continue;
        }




        // --- not in heap yet --- 
        if (k < X){
            heap_push(heap, &k, (Hit){ d_eq, i });
            if (k == X){
                cutoff_km      = heap_top(heap).dist_km * (1.0 + EPSILON);
                lat_thresh_deg = cutoff_km / 111.0;
                lon_thresh_deg = use_lon ? cutoff_km / (111.0 * coslat) : 360.0;
            }
        } else {
            // heap full: candidate must be competitive vs current worst 
            double worst = heap_top(heap).dist_km;
            if (d_eq <= worst * (1.0 + EPSILON)){
                heap_replace_top(heap, k, (Hit){ d_eq, i });
                cutoff_km      = heap_top(heap).dist_km * (1.0 + EPSILON);
                lat_thresh_deg = cutoff_km / 111.0;
                lon_thresh_deg = use_lon ? cutoff_km / (111.0 * coslat) : 360.0;
            }
        }
    }

    if (k == 0){
        puts("No candidates after streaming.");
        free(heap);
        return;
    }

    // refine: just reuse equirectangular distances and sort 
    Hit *refined = (Hit*)malloc((size_t)k * sizeof(Hit));
    if (!refined){
        fprintf(stderr,"OOM phase4 refined\n");
        free(heap);
        return;
    }
    for (int j=0; j<k; ++j){
        int idx = heap[j].idx;
        refined[j].idx     = idx;
        refined[j].dist_km = equirectangular(ref_lat, ref_lon, arr[idx].lat, arr[idx].lon);
    }
    qsort(refined, (size_t)k, sizeof(Hit), cmp_hit_dist);

    int *idx_refined = (int*)malloc((size_t)k * sizeof(int));
    if (!idx_refined){
        fprintf(stderr,"OOM phase4 idx_refined\n");
        free(refined);
        free(heap);
        return;
    }
    for (int i=0; i<k; ++i) {
        int idx = refined[i].idx;
        idx_refined[i] = idx;
        arr[idx].metric = refined[i].dist_km; // store equirectangular km 
    }

    print_top_from_idx(arr, idx_refined, k, X, 4, ref_lat, ref_lon,
                       "Phase 4 (Adaptive Pruning, Equirectangular)");
    write_csv_from_idx("output/phase4-output-serial.csv", arr, idx_refined, k, X, 4, ref_lat, ref_lon);

    free(idx_refined);
    free(refined);
    free(heap);

    double t1 = wall_time();
    double elapsed_ms = (t1 - t0) * 1000.0;
    printf("Bounding-box survivors: %d of %d\n", bbox_survivors, count);
    printf("Post-equirectangular survivors: %d\n", k);
    printf("Phase 4 runtime: %.3f ms\n", elapsed_ms);
    log_timing(LOG_FILE,"Phase4_PruneRefine",count,1,elapsed_ms);
}

// ---------- main ---------- 
int main(int argc,char *argv[]){
    Args args = parse_args(argc, argv);

    if (args.phase < 1 || args.phase > 4){
        fprintf(stderr,
            "Usage: %s --phase {1|2|3|4} [lat lon X] [--data-file <path>] [--use-rows N]\n", argv[0]);
        fprintf(stderr,
            "Example: %s --phase=3 30.4733 -87.1866 25 --data-file input/myfile.csv --use-rows 100000 \n", argv[0]);
        return 1;
    }

    const char *auto_name =
        (args.phase==1) ? "output/phase1-output-serial.csv" :
        (args.phase==2) ? "output/phase2-output-serial.csv" :
        (args.phase==3) ? "output/phase3-output-serial.csv" :
                          "output/phase4-output-serial.csv";

    if (ensure_dir("output") != 0){
        perror("ensure_dir(\"output\")");
    }

    printf("Phase: %d\n", args.phase);
    printf("Output: %s\n", auto_name);
    printf("Ref: (%.4f, %.4f), TopX=%d\n", args.ref_lat, args.ref_lon, args.topX);
    printf("Data file: %s\n", args.data_file);
    if (args.use_rows > 0) {
        printf("Use rows (max): %d\n", args.use_rows);
    }

    int capacity = 0;
    int count = 0;
    Aircraft *aircrafts = NULL;

    int rc;
    // If data_file contains *, ?, or [ ], treat it as a glob pattern
    
    if (strpbrk(args.data_file, "*?[]") != NULL) {
        rc = load_csv_pattern(args.data_file, &aircrafts, &count, &capacity);
    } else {
        rc = load_csv(args.data_file, &aircrafts, &count, &capacity);
    }

    if (rc != 0) {
        fprintf(stderr, "Failed to load from '%s'\n", args.data_file);
        if (aircrafts) free(aircrafts);
        return 1;
    }

    printf("Loaded %d records\n", count);

    if (args.use_rows > 0 && count > args.use_rows) {
        printf("Capping to first %d rows due to --use-rows\n", args.use_rows);
        count = args.use_rows;
    }

    if (count == 0) {
        free(aircrafts);
        return 0;
    }  

    if      (args.phase==1) phase1(aircrafts,count,args.topX,args.ref_lat,args.ref_lon);
    else if (args.phase==2) phase2(aircrafts,count,args.topX,args.ref_lat,args.ref_lon);
    else if (args.phase==3) phase3(aircrafts,count,args.topX,args.ref_lat,args.ref_lon);
    else                    phase4(aircrafts,count,args.topX,args.ref_lat,args.ref_lon);

    free(aircrafts);
    return 0;
}

