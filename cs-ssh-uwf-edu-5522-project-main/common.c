#define _POSIX_C_SOURCE 200809L
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "common.h"
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <time.h>
#include <stdio.h>
#include <glob.h>
#include <limits.h>

/* These are implemented in the drivers; declare so we can call them here. */
double full_haversine(double lat1,double lon1,double lat2,double lon2);
double distance_from_half(double a);


/* -------- timing / fs ---------- */
double wall_time(void){
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int ensure_dir(const char *path){
    struct stat st;
    if (stat(path, &st) == 0) {
        if (S_ISDIR(st.st_mode)) return 0;
        errno = ENOTDIR;
        return -1;
    }
    if (mkdir(path, 0755) == 0) return 0;
    if (errno == EEXIST) return 0;
    return -1;
}

/* -------- strings ---------- */
void rstrip_newline(char *s){
    if(!s) return;
    size_t n = strlen(s);
    while(n && (s[n-1]=='\n' || s[n-1]=='\r')) s[--n] = '\0';
}

void rtrim_spaces(char *s){
    if(!s) return;
    size_t n = strlen(s);
    while(n && isspace((unsigned char)s[n-1])) s[--n] = '\0';
}

int split_csv_simple(char *line, char *fields[], int max_fields){
    int n = 0;
    char *p = line, *start = line;
    while(*p && n < max_fields){
        if(*p == ','){
            *p = '\0';
            fields[n++] = start;
            start = p + 1;
        }
        p++;
    }
    if(n < max_fields) fields[n++] = start;
    return n;
}

bool parse_double_strict(const char *s, double *out){
    if(!s || !*s) return false;
    errno = 0;
    char *end = NULL;
    double v = strtod(s, &end);
    if(errno || end == s || *end != '\0' || !isfinite(v)) return false;
    *out = v;
    return true;
}

double parse_double_or0(const char *s){
    if(!s || !*s) return 0.0;
    char *end = NULL;
    double v = strtod(s, &end);
    return (end && *end == '\0') ? v : 0.0;
}

const char* dash_if_empty(const char *s){
    return (s && *s) ? s : "-";
}



/* Helper: from a sorted idx[] (by metric), pick top X unique aircraft (by ICAO24).
 * Returns the number of selected aircraft (<= X), stored in out_idx[].
 */
static int select_top_unique_icao24(const Aircraft R[],
                                    const int *idx_sorted,
                                    int N, int X,
                                    int *out_idx)
{
    if (N <= 0 || X <= 0) return 0;
    int selected = 0;

    for (int i = 0; i < N && selected < X; ++i) {
        int k = idx_sorted[i];
        const char *icao = R[k].icao24;
        bool already = false;

        /* Check if we've already picked this ICAO */
        for (int j = 0; j < selected; ++j) {
            int k2 = out_idx[j];
            if (strcmp(R[k2].icao24, icao) == 0) {
                already = true;
                break;
            }
        }
        if (!already) {
            out_idx[selected++] = k;
        }
    }
    return selected;
}




/* -------- CSV load / print / logs ---------- */

int load_csv_pattern(const char *pattern,
                     Aircraft **arr, int *count, int *capacity)
{
    glob_t g;
    memset(&g, 0, sizeof(g));

    int rv = glob(pattern, 0, NULL, &g);
    if (rv != 0) {
        if (rv == GLOB_NOMATCH) {
            fprintf(stderr,
                    "No files matched pattern '%s'\n", pattern);
        } else {
            fprintf(stderr,
                    "glob() failed for pattern '%s' (code %d)\n",
                    pattern, rv);
        }
        globfree(&g);
        return -1;
    }

    // glob() returns file names in sorted order by default.
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        const char *fname = g.gl_pathv[i];
        fprintf(stderr, "Loading %s\n", fname);
        int rc = load_csv(fname, arr, count, capacity);
        if (rc != 0) {
            fprintf(stderr, "Failed to load %s (rc=%d)\n", fname, rc);
            globfree(&g);
            return rc;
        }
    }

    globfree(&g);
    return 0;
}


int load_csv(const char *filename, Aircraft **arr, int *count, int *capacity){
    FILE *fp = NULL;
    int is_pipe = 0;

    /* Detect .gz extension (simple suffix check) */
    size_t len = strlen(filename);
    if (len >= 3 && strcmp(filename + len - 3, ".gz") == 0) {
        /* Read via gzip -cd so we still get a FILE* stream */
        char cmd[PATH_MAX + 32];
        /* Quote with single quotes to handle spaces in path */
        snprintf(cmd, sizeof(cmd), "gzip -cd -- '%s'", filename);

        fp = popen(cmd, "r");
        if (!fp) {
            perror("popen gzip");
            return -1;
        }
        is_pipe = 1;
    } else {
        fp = fopen(filename,"r");
        if (!fp) {
            perror("CSV open");
            return -1;
        }
    }

    if ((*arr == NULL) || (*capacity <= 0)) {
        *capacity = (*capacity > 0) ? *capacity : INIT_CAPACITY;
        *arr = (Aircraft*)malloc((size_t)(*capacity) * sizeof(Aircraft));
        if(!*arr){
            perror("malloc");
            if (is_pipe) pclose(fp); else fclose(fp);
            return -1;
        }
        *count = 0;
    }

    char line[MAX_LINE];

    int idx_time=-1, idx_icao=-1, idx_lat=-1, idx_lon=-1,
        idx_vel=-1, idx_head=-1, idx_callsign=-1, idx_baro=-1;

    /* Read header line */
    if (!fgets(line, MAX_LINE, fp)){
        if (is_pipe) pclose(fp); else fclose(fp);
        return 0;
    }
    rstrip_newline(line);

    {
        char *fields[64];
        int nf = split_csv_simple(line, fields, 64);

        /* Skip our own output files (rank header) */
        if (nf > 0 && !strcmp(fields[0], "rank")) {
            if (is_pipe) pclose(fp); else fclose(fp);
            return 0;
        }

        for (int i=0;i<nf;i++){
            if      (!strcmp(fields[i],"time"))          idx_time = i;
            else if (!strcmp(fields[i],"icao24"))        idx_icao = i;
            else if (!strcmp(fields[i],"lat"))           idx_lat  = i;
            else if (!strcmp(fields[i],"lon"))           idx_lon  = i;
            else if (!strcmp(fields[i],"velocity"))      idx_vel  = i;
            else if (!strcmp(fields[i],"heading"))       idx_head = i;
            else if (!strcmp(fields[i],"callsign"))      idx_callsign = i;
            else if (!strcmp(fields[i],"baroaltitude"))  idx_baro = i;
        }

        /* Fallback for legacy 8-col input (positional) */
        if (idx_time<0 && nf>=8){
            idx_time=0; idx_icao=1; idx_lat=2; idx_lon=3;
            idx_vel=4; idx_head=5; idx_callsign=6; idx_baro=7;
        }
    }

    while (fgets(line, MAX_LINE, fp)){
        rstrip_newline(line);

        if(*count >= *capacity){
            int newcap = (*capacity) * 2;
            if (newcap < 1) newcap = INIT_CAPACITY;
            Aircraft *tmp = (Aircraft*)realloc(*arr, (size_t)newcap * sizeof(Aircraft));
            if(!tmp){
                perror("realloc");
                if (is_pipe) pclose(fp); else fclose(fp);
                return -1;
            }
            *arr = tmp;
            *capacity = newcap;
        }

        char *fields[64];
        int nf = split_csv_simple(line, fields, 64);

        Aircraft ac; memset(&ac, 0, sizeof(ac));

        if (idx_time>=0 && idx_time<nf){
            strncpy(ac.time, fields[idx_time], 31);
            ac.time[31]='\0';
            rtrim_spaces(ac.time);
        }
        if (idx_icao>=0 && idx_icao<nf){
            strncpy(ac.icao24, fields[idx_icao], 15);
            ac.icao24[15]='\0';
            rtrim_spaces(ac.icao24);
        }
        if (idx_callsign>=0 && idx_callsign<nf){
            strncpy(ac.callsign, fields[idx_callsign], 15);
            ac.callsign[15]='\0';
            rtrim_spaces(ac.callsign);
        }
        if (idx_vel>=0    && idx_vel<nf)   ac.velocity     = parse_double_or0(fields[idx_vel]);
        if (idx_head>=0   && idx_head<nf)  ac.heading      = parse_double_or0(fields[idx_head]);
        if (idx_baro>=0   && idx_baro<nf)  ac.baroaltitude = parse_double_or0(fields[idx_baro]);

        double lat = 0.0, lon = 0.0;
        bool lat_ok = (idx_lat>=0 && idx_lat<nf) && parse_double_strict(fields[idx_lat], &lat);
        bool lon_ok = (idx_lon>=0 && idx_lon<nf) && parse_double_strict(fields[idx_lon], &lon);
        if (!lat_ok || !lon_ok || lat < -90.0 || lat > 90.0 || lon < -180.0 || lon > 180.0){
            continue;  /* ignore this record */
        }
        ac.lat = lat; ac.lon = lon;

        (*arr)[(*count)++] = ac;
    }

    if (is_pipe) pclose(fp); else fclose(fp);
    return 0;
}


void print_top_from_idx(const Aircraft R[], const int *idx, int N, int X,
                        int phase, double ref_lat, double ref_lon, const char *label){

    if (N <= 0 || X <= 0) {
        printf("\n--- %s ---\n(no results)\n", label ? label : "");
        return;
    }

    /* First, select up to X unique aircraft by ICAO24 from the sorted idx[] */
    int *sel = (int*)malloc(sizeof(int) * (size_t)X);
    if (!sel) {
        fprintf(stderr, "OOM in print_top_from_idx\n");
        return;
    }
    int used = select_top_unique_icao24(R, idx, N, X, sel);
    if (used <= 0) {
        printf("\n--- %s ---\n(no unique aircraft)\n", label ? label : "");
        free(sel);
        return;
    }
    if (used < X) X = used;  /* clamp */

    printf("\n--- %s ---\n", label ? label : "Top");
    printf("  %-4s %-8s %-8s %9s %10s %8s %10s\n",
           "Rank","ICAO","Callsign","Lat(deg)","Lon(deg)","Alt","Dist(km)");
    printf("  %-4s %-8s %-10s %9s %10s %8s %10s\n",
           "----","--------","--------","---------","----------","--------","----------");

    for (int i = 0; i < X; i++){
        int k = sel[i];
        double dist = 0.0;

        if      (phase == 1) {
            /* Phase 1: metric holds Vincenty distance (km) */
            dist = R[k].metric;
        } else if (phase == 2) {
            /* Phase 2: metric holds full-haversine distance (km) */
            dist = R[k].metric;
        } else if (phase == 3) {
            /* Phase 3: metric holds half-haversine "a" -> convert to km */
            dist = distance_from_half(R[k].metric);
        } else if (phase == 4) {
            /* Phase 4: metric holds equirectangular distance (km) */
            dist = R[k].metric;
        } else {
            /* Fallback: compute full haversine directly */
            dist = full_haversine(ref_lat, ref_lon, R[k].lat, R[k].lon);
        }

        printf("  %4d %-8s %-8s %9.4f %10.4f %8.0f %10.4f\n",
               i+1, dash_if_empty(R[k].icao24), dash_if_empty(R[k].callsign),
               R[k].lat, R[k].lon, R[k].baroaltitude, dist);
    }

    free(sel);
}




void write_csv_from_idx(const char *filename, const Aircraft R[], const int *idx,
                        int N, int X, int phase, double ref_lat, double ref_lon)
{
    if (N <= 0 || X <= 0) {
        printf("No rows to write.\n");
        return;
    }

    int *sel = (int*)malloc(sizeof(int) * (size_t)X);
    if (!sel) {
        fprintf(stderr, "OOM in write_csv_from_idx\n");
        return;
    }

    int used = select_top_unique_icao24(R, idx, N, X, sel);
    if (used <= 0) {
        printf("No unique aircraft to write.\n");
        free(sel);
        return;
    }
    if (used < X) X = used;  /* clamp */

    FILE *fp = fopen(filename,"w");
    if(!fp){
        perror("CSV write");
        free(sel);
        return;
    }
    fprintf(fp,"rank,time,icao24,callsign,lat,lon,baroaltitude,dist_km\n");

    for(int i=0;i<X;i++){
        int k = sel[i];
        double dist = 0.0;
        if      (phase == 1) {
            /* Phase 1: metric holds Vincenty distance in km */
            dist = R[k].metric;
        } else if (phase == 2) {
            /* Phase 2: metric holds full-haversine distance in km */
            dist = R[k].metric;
        } else if (phase == 3) {
            /* Phase 3: metric holds half-haversine "a" -> convert to km */
            dist = distance_from_half(R[k].metric);
        } else if (phase == 4) {
            /* Phase 4: metric holds equirectangular distance in km */
            dist = R[k].metric;
        } else {
            dist = full_haversine(ref_lat, ref_lon, R[k].lat, R[k].lon);
        }
        fprintf(fp,"%d,%s,%s,%s,%.6f,%.6f,%.1f,%.3f\n",
                i+1,R[k].time,R[k].icao24,R[k].callsign,
                R[k].lat,R[k].lon,R[k].baroaltitude,dist);
    }

    fclose(fp);
    free(sel);
    printf("Results written to %s\n",filename);
}



void log_timing(const char *filename,const char *phase_name,
                int count,int cores,double runtime_ms){
    FILE *fp=fopen(filename,"a");
    if(!fp){ perror("Timing log"); return; }
    fprintf(fp,"%s,%d,%d,%.6f\n",phase_name,count,cores,runtime_ms);
    fclose(fp);
}

/* -------- heap for phase3/4 ---------- */
void heap_sift_up(Hit h[], int i){
    while(i>0){
        int p=(i-1)/2;
        if(h[p].dist_km>=h[i].dist_km) break;
        Hit t=h[p]; h[p]=h[i]; h[i]=t; i=p;
    }
}

void heap_sift_down(Hit h[], int n, int i){
    for(;;){
        int l=2*i+1,r=l+1,m=i;
        if(l<n&&h[l].dist_km>h[m].dist_km) m=l;
        if(r<n&&h[r].dist_km>h[m].dist_km) m=r;
        if(m==i) break;
        Hit t=h[i]; h[i]=h[m]; h[m]=t; i=m;
    }
}

void heap_push(Hit h[], int *n, Hit v){ h[(*n)++]=v; heap_sift_up(h,*n-1); }
Hit  heap_top(Hit h[]){ return h[0]; }
void heap_replace_top(Hit h[], int n, Hit v){ (void)n; h[0]=v; heap_sift_down(h,n,0); }

/* ---------- CLI ---------- */
Args parse_args(int argc, char *argv[]){
    Args a; memset(&a,0,sizeof(a));
    a.phase     = 0;
    a.output    = NULL; /* unused */
    a.ref_lat   = DEFAULT_LAT;
    a.ref_lon   = DEFAULT_LON;
    a.topX      = DEFAULT_X;
    a.data_file = DATA_FILE;
    a.use_rows  = 0;    /* 0 = no limit */

    for (int i=1; i<argc; ++i){
        if (starts_with(argv[i],"--phase=")) {
            a.phase = atoi(argv[i]+8);
        } else if (!strcmp(argv[i],"--phase") && i+1<argc) {
            a.phase = atoi(argv[++i]);
        } else if (starts_with(argv[i],"--data-file=")) {
            a.data_file = argv[i] + 12;
        } else if (!strcmp(argv[i],"--data-file") && i+1<argc) {
            a.data_file = argv[++i];

        } else if (starts_with(argv[i],"--use-rows=")) {
            a.use_rows = (int)atol(argv[i] + 11);
        } else if (!strcmp(argv[i],"--use-rows") && i+1<argc) {
            a.use_rows = (int)atol(argv[++i]);
        }
    }

    /* Trailing [lat lon X] */
    if (argc>=4){
        char *e1=NULL,*e2=NULL,*e3=NULL;
        double tlat = strtod(argv[argc-3], &e1);
        double tlon = strtod(argv[argc-2], &e2);
        long   tX   = strtol(argv[argc-1], &e3, 10);
        if (e1 && *e1=='\0' && e2 && *e2=='\0' && e3 && *e3=='\0'){
            a.ref_lat = tlat; a.ref_lon = tlon; a.topX = (int)tX;
        }
    }
    return a;
}



/* Robust Vincenty inverse; falls back to haversine if it fails to converge. (near poles)*/
double vincenty_distance_km(double lat1, double lon1,
                                   double lat2, double lon2)
{
    if (lat1 == lat2 && lon1 == lon2) return 0.0;

    double phi1 = deg2rad(lat1);
    double phi2 = deg2rad(lat2);
    double U1 = atan((1.0 - WGS84_F) * tan(phi1));
    double U2 = atan((1.0 - WGS84_F) * tan(phi2));
    double sinU1 = sin(U1), cosU1 = cos(U1);
    double sinU2 = sin(U2), cosU2 = cos(U2);

    double L = deg2rad(lon_diff(lon2, lon1));
    double lambda = L, lambda_prev;
    double sinSigma, cosSigma, sigma;
    double sinAlpha, cos2Alpha, cos2SigmaM, C;

    for (int iter = 0; iter < 100; ++iter) {
        double sinLambda = sin(lambda);
        double cosLambda = cos(lambda);

        double t1 = cosU2 * sinLambda;
        double t2 = cosU1 * sinU2 - sinU1 * cosU2 * cosLambda;
        sinSigma  = sqrt(t1*t1 + t2*t2);
        if (sinSigma == 0.0) return 0.0;

        cosSigma  = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
        sigma     = atan2(sinSigma, cosSigma);

        sinAlpha  = (cosU1 * cosU2 * sinLambda) / sinSigma;
        cos2Alpha = 1.0 - sinAlpha * sinAlpha;

        if (cos2Alpha != 0.0)
            cos2SigmaM = cosSigma - 2.0 * sinU1 * sinU2 / cos2Alpha;
        else
            cos2SigmaM = 0.0; /* equatorial */

        C = (WGS84_F / 16.0) * cos2Alpha * (4.0 + WGS84_F * (4.0 - 3.0 * cos2Alpha));

        lambda_prev = lambda;
        lambda = L + (1.0 - C) * WGS84_F * sinAlpha *
                 (sigma + C * sinSigma *
                  (cos2SigmaM + C * cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM)));

        if (fabs(lambda - lambda_prev) < 1e-12) break;

        if (iter == 99) {
            /* fallback: spherical distance */
            return full_haversine(lat1, lon1, lat2, lon2);
        }
    }

    double uSq = cos2Alpha * ((WGS84_A*WGS84_A - WGS84_B*WGS84_B) / (WGS84_B*WGS84_B));
    double A = 1.0 + uSq/16384.0 * (4096.0 + uSq * (-768.0 + uSq * (320.0 - 175.0*uSq)));
    double B =       uSq/1024.0  * ( 256.0 + uSq * (-128.0 + uSq * ( 74.0 -  47.0*uSq)));
    double deltaSigma = B * sinSigma *
        (cos2SigmaM + (B/4.0) * (cosSigma * (-1.0 + 2.0 * cos2SigmaM * cos2SigmaM) -
         (B/6.0) * cos2SigmaM * (-3.0 + 4.0 * sinSigma * sinSigma)
                            * (-3.0 + 4.0 * cos2SigmaM * cos2SigmaM)));

    double s_m = WGS84_B * A * (sigma - deltaSigma); /* meters */
    return s_m * 1e-3; /* km */
}

