/* C caller example: replay RTCM3 binary or NMEA text; output final snapshot.
   od_replay rtcm|nmea INPUT [observation_capacity] */
#include "orbit_determination.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    od_context_t *ctx = NULL;
    od_config_t config = od_default_config();
    od_j2000_state_t state;
    uint8_t buffer[257]; /* Deliberately unrelated to frame boundaries. */
    size_t length;
    uint64_t accepted = 0, rejected = 0;
    FILE *input;
    od_status_t status;
    int is_rtcm;
    if (argc < 3 || argc > 4 || (strcmp(argv[1], "rtcm") && strcmp(argv[1], "nmea"))) {
        fprintf(stderr, "usage: %s rtcm|nmea INPUT [capacity]\n", argv[0]);
        return 2;
    }
    is_rtcm = strcmp(argv[1], "rtcm") == 0;
    if (argc == 4) {
        char *end;
        unsigned long capacity = strtoul(argv[3], &end, 10);
        if (!argv[3][0] || *end || capacity < 3 || capacity > 65536) return 2;
        config.observation_capacity = (uint32_t)capacity;
    }
    status = od_create(&config, &ctx);
    if (status != OD_OK) { fprintf(stderr, "%s\n", od_status_string(status)); return 1; }
    input = fopen(argv[2], "rb");
    if (!input) { perror(argv[2]); od_destroy(ctx); return 1; }
    while ((length = fread(buffer, 1, sizeof(buffer), input)) != 0) {
        od_feed_info_t info;
        status = is_rtcm ? od_feed_rtcm(ctx, buffer, length, &info)
                         : od_feed_nmea(ctx, (const char *)buffer, length, &info);
        accepted += info.accepted_observations;
        rejected += info.rejected_records;
        if (status != OD_OK) break;
    }
    if (ferror(input)) { perror("read"); status = OD_ERROR_INTERNAL; }
    fclose(input);
    fprintf(stderr, "accepted=%" PRIu64 ", rejected=%" PRIu64 ", capacity=%u\n",
            accepted, rejected, config.observation_capacity);
    if (status == OD_OK) status = od_get_orbit(ctx, &state);
    od_destroy(ctx);
    if (status != OD_OK) { fprintf(stderr, "%s\n", od_status_string(status)); return 1; }
    puts("timestamp_utc_ms,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps");
    printf("%" PRId64 ",%.9f,%.9f,%.9f,%.9f,%.9f,%.9f\n", state.timestamp_utc_ms,
        state.position_m[0], state.position_m[1], state.position_m[2],
        state.velocity_mps[0], state.velocity_mps[1], state.velocity_mps[2]);
    return 0;
}
