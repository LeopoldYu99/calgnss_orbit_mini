#ifndef CALGNSS_ORBIT_MINI_H
#define CALGNSS_ORBIT_MINI_H

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
} cg_observation_t;

typedef struct cg_options_t {
    int degree;
    /* <= 0 means no fixed future extrapolation limit. */
    double max_extrapolation_seconds;
    /* <= 0 means use all cached history for future extrapolation initial-state fitting. */
    double extrapolation_history_seconds;
    double propagation_step_seconds;
    int enable_orbit_phase_correction;
} cg_options_t;

#define CG_MAX_DEGREE 16
#define CG_DEFAULT_OBSERVATION_CAPACITY 10

typedef struct cg_context_t cg_context_t;

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

int cg_precompute_observations(cg_observation_t *observations, size_t count);

int cg_context_create(
    cg_context_t **out_context,
    cg_observation_t *observation_buffer,
    size_t capacity,
    const cg_options_t *options);

void cg_context_reset(cg_context_t *context);

void cg_context_destroy(cg_context_t *context);

size_t cg_context_count(const cg_context_t *context);

size_t cg_context_capacity(const cg_context_t *context);

int cg_context_push(
    cg_context_t *context,
    const cg_observation_t *observation);

int cg_context_query_state(
    cg_context_t *context,
    const cg_time_t *query_time_utc,
    cg_state_t *out_state);

const char *cg_status_string(int status);

#ifdef __cplusplus
}
#endif

#endif
