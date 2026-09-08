#include "orbit_determination.h"

int main(void)
{
    od_context_t *ctx = 0;
    od_feed_info_t info;
    od_j2000_state_t state;
    static const char nmea[] =
        "$GNRMC,110000.00,A,,,,,,,150525,,,A\r\n"
        "$GNGGA,110000.00,0012.6157446,S,16642.1045030,E,1,12,0.8,500000.0,M,0.0,M,,\r\n"
        "$GNGGA,110002.00,0012.5495000,S,16642.0900000,E,1,12,0.8,500000.0,M,0.0,M,,\r\n"
        "$GNGGA,110004.00,0012.4833000,S,16642.0755000,E,1,12,0.8,500000.0,M,0.0,M,,\r\n"
        "$GNGGA,110006.00,0012.4171000,S,16642.0610000,E,1,12,0.8,500000.0,M,0.0,M,,\r\n";
    if (od_create(0, &ctx) != OD_OK) return 1;
    if (od_feed_nmea(ctx, nmea, sizeof(nmea)-1, &info) != OD_OK) return 2;
    if (info.accepted_observations != 4 || od_get_orbit(ctx, &state) != OD_OK) return 3;
    od_destroy(ctx);
    return 0;
}
