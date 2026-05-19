#ifndef CALGNSS_H
#define CALGNSS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct cg_time_t {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    double jd_utc;
    double unix_seconds;
} cg_time_t;

typedef struct cg_vec3_t {
    double x;
    double y;
    double z;
} cg_vec3_t;

typedef struct cg_state_t {
    cg_time_t time_utc;
    cg_vec3_t r_j2000_m;
    cg_vec3_t v_j2000_mps;
} cg_state_t;

typedef struct cg_observation_t {
    cg_time_t time_utc;
    double lat_deg;
    double lon_deg;
    double alt_m;
    int has_rates;
    double lat_rate_degps;
    double lon_rate_degps;
    double alt_rate_mps;
    cg_vec3_t r_ecef_m;
    cg_vec3_t v_ecef_mps;
    cg_vec3_t r_j2000_m;
    cg_vec3_t v_j2000_mps;
} cg_observation_t;

typedef struct cg_options_t {
    int degree;
    double fit_window_minutes;
    /* <= 0 means no fixed future extrapolation limit. */
    double max_extrapolation_seconds;
    double propagation_step_seconds;
    int enable_orbit_phase_correction;
} cg_options_t;

typedef enum cg_status_t {
    CG_OK = 0,
    CG_ERR_INVALID_ARGUMENT = -1,
    CG_ERR_PARSE = -2,
    CG_ERR_IO = -3,
    CG_ERR_NO_MEMORY = -4,
    CG_ERR_RANGE = -5,
    CG_ERR_FIT = -6
} cg_status_t;

cg_options_t cg_default_options(void);

int cg_parse_time(const char *text, cg_time_t *out_time);
void cg_format_time_iso(const cg_time_t *time_utc, char *buffer, size_t buffer_size);
double cg_seconds_between(const cg_time_t *a, const cg_time_t *b);
cg_time_t cg_time_add_seconds(const cg_time_t *time_utc, double seconds);

int cg_load_lla_csv(const char *path, cg_observation_t **out_observations, size_t *out_count);
void cg_free_observations(cg_observation_t *observations);
int cg_precompute_observations(cg_observation_t *observations, size_t count);

int cg_query_state(
    const cg_observation_t *observations,
    size_t count,
    const cg_time_t *query_time_utc,
    const cg_options_t *options,
    cg_state_t *out_state);

int cg_export_j2000_csv(
    const char *input_lla_csv,
    const char *output_j2000_csv,
    const cg_time_t *start_time_utc,
    const cg_time_t *end_time_utc,
    double step_seconds,
    const cg_options_t *options);

const char *cg_status_string(int status);

#ifdef __cplusplus
}
#endif

#endif
