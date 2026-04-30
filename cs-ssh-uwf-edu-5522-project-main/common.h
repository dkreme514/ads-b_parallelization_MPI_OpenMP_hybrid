#ifndef COMMON_H
#define COMMON_H

#include <stdbool.h>
#include <stddef.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef EARTH_RADIUS_KM
#define EARTH_RADIUS_KM 6371.0
#endif


#ifndef WGS84_A
#define WGS84_A 6378137.0                 // meters 
#endif

#ifndef WGS84_F
#define WGS84_F (1.0/298.257223563)
#endif

#ifndef WGS84_B
#define WGS84_B (WGS84_A * (1.0 - WGS84_F))
#endif

#ifndef MAX_LINE
#define MAX_LINE 2048       // Max size of line fetched from CSV 
#endif

#ifndef INIT_CAPACITY
#define INIT_CAPACITY 512   // Initial Aircraft array size.. doubles as it fills 
#endif

// CLI / app defaults (shared) 
#ifndef DEFAULT_LAT
#define DEFAULT_LAT 30.5468
#endif
#ifndef DEFAULT_LON
#define DEFAULT_LON -87.2174
#endif
#ifndef DEFAULT_X
#define DEFAULT_X 25
#endif
#ifndef DATA_FILE
#define DATA_FILE "input/states_2022-*.csv"
#endif
#define LOG_FILE "input/timings.csv"

// ---------- types shared across units ---------- 
typedef struct {
    char time[32];
    char icao24[16];
    double lat, lon;
    double velocity;
    double heading;
    double baroaltitude;
    char callsign[16];
    double metric;   // phase-dependent metric used for sorting 
    double dist_km;  // final computed distance for printing/CSV 
} Aircraft;

typedef struct { double dist_km; int idx; } Hit;

typedef struct {
    int   phase;         // 1,2,3
    char *output;        // (ignored; output is phase-specific)
    double ref_lat;      // optional override
    double ref_lon;
    int    topX;         // optional override
    const char *data_file; // CLI overridable input path
    int    use_rows;     // 0 = no limit, >0 = cap number of rows used
} Args;


// ---------- tiny helpers ---------- 
static inline int starts_with(const char *s, const char *p) {
    if (!s || !p) return 0;
    while (*p) {
        if (*s++ != *p++) return 0;
    }
    return 1;
}

static inline double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

static inline double lon_diff(double lon1, double lon2) {
    double dlon = lon1 - lon2;
    while (dlon > 180.0)  dlon -= 360.0;
    while (dlon < -180.0) dlon += 360.0;
    return dlon;
}



// distance helpers implemented in the drivers 
double full_haversine(double lat1,double lon1,double lat2,double lon2);
double half_haversine(double lat1,double lon1,double lat2,double lon2);
double distance_from_half(double a);
double vincenty_distance_km(double lat1, double lon1,double lat2, double lon2);

// ---------- utilities provided by common.c ---------- 
double wall_time(void);
int    ensure_dir(const char *path);

void rstrip_newline(char *s);
void rtrim_spaces(char *s);
int  split_csv_simple(char *line, char *fields[], int max_fields);
bool parse_double_strict(const char *s, double *out);
double parse_double_or0(const char *s);
const char* dash_if_empty(const char *s);

// CSV I/O, printing, logging (these will call distance helpers provided in drivers) 
int load_csv_pattern(const char *pattern, Aircraft **arr, int *count, int *capacity);
int  load_csv(const char *filename, Aircraft **arr, int *count, int *capacity);
void print_top_from_idx(const Aircraft R[], const int *idx, int N, int X,
                        int phase, double ref_lat, double ref_lon, const char *label);
void write_csv_from_idx(const char *filename, const Aircraft R[], const int *idx,
                        int N, int X, int phase, double ref_lat, double ref_lon);
void log_timing(const char *filename,const char *phase_name,
                int count,int cores,double runtime);

// compact max-heap used by phase 3 (provided in common.c) 
void heap_sift_up(Hit h[], int i);
void heap_sift_down(Hit h[], int n, int i);
void heap_push(Hit h[], int *n, Hit v);
Hit  heap_top(Hit h[]);
void heap_replace_top(Hit h[], int n, Hit v);

// CLI parser shared by both drivers 
Args parse_args(int argc, char *argv[]);

// NOTE: No declarations here for deg2rad/lon_diff/haversine/equirectangular
// and no declarations for sorting. Those live inside each driver. 

#endif // COMMON_H 

