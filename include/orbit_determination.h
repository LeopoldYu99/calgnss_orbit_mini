#ifndef ORBIT_DETERMINATION_H
#define ORBIT_DETERMINATION_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct od_context od_context_t;
typedef int32_t od_status_t;

enum {
    OD_OK = 0,
    OD_NOT_READY = 1,
    OD_ERROR_ARGUMENT = -1,
    OD_ERROR_MEMORY = -2,
    OD_ERROR_SOURCE = -3,
    OD_ERROR_FIT = -4,
    OD_ERROR_INTERNAL = -5
};

typedef struct od_config {
    /* Ring capacity in position observations, NOT bytes. Range: 2..65536.
       Once full, each new observation replaces the oldest. Default: 30. */
    uint32_t observation_capacity;
    /* Polynomial degree: 1..16, strictly less than capacity. Default: 3.
       Output becomes ready after fit_degree + 1 observations. */
    uint32_t fit_degree;
    /* UTC Unix milliseconds for resolving GGA time without RMC.
       -1: no reference; require RMC. Reference must be within 12 hours of
       the first GGA. Refresh with RMC after gaps >= 12 hours. */
    int64_t nmea_reference_utc_ms;
    /* Added to RTCM 1019's 10-bit GPS week. Multiple of 1024, <= 8192.
       Default 2048: GPS weeks 2048..3071. Never reads host clock. */
    uint32_t gps_week_rollover;
} od_config_t;

typedef struct od_j2000_state {
    int64_t timestamp_utc_ms;
    double position_m[3];
    double velocity_mps[3];
} od_j2000_state_t;

typedef struct od_feed_info {
    uint32_t accepted_observations;
    uint32_t rejected_records;
    uint32_t buffered_observations;
} od_feed_info_t;

od_config_t od_default_config(void);
/* config == NULL selects defaults. On failure *out_ctx is NULL. */
od_status_t od_create(const od_config_t *config, od_context_t **out_ctx);
void od_destroy(od_context_t *ctx);
/* Clear observations, decoders, ephemerides and input-source selection;
   retain configuration. */
od_status_t od_reset(od_context_t *ctx);

/* Two streaming inputs. Supports split/concatenated records. Max chunk:
   1 MiB. length == 0 is a no-op. out_info may be NULL. OD_OK means input
   processed/buffered; inspect counts to see accepted/rejected records.
   Frames awaiting ephemerides and valid unsupported messages add no point.
   Malformed, missing-date and non-increasing observations are discarded.
   One source per context until reset. Caller serializes calls per context. */
od_status_t od_feed_rtcm(od_context_t *ctx, const uint8_t *data,
                         size_t length, od_feed_info_t *out_info);
od_status_t od_feed_nmea(od_context_t *ctx, const char *data,
                         size_t length, od_feed_info_t *out_info);

/* Latest observation epoch only: no query time or future extrapolation.
   Fit uses all currently buffered positions. Output is unchanged on error
   or OD_NOT_READY. Repeated calls without new observations return the same
   snapshot. Feed/get per epoch if intermediate snapshots are required. */
od_status_t od_get_orbit(od_context_t *ctx, od_j2000_state_t *out_state);
const char *od_status_string(od_status_t status);

#ifdef __cplusplus
}
#endif
#endif
