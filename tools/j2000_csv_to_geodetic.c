/* Convert captured J2000 orbit points back to ITRS/ECEF and WGS-84 geodetic.
 *
 * This diagnostic intentionally includes the compact orbit implementation so
 * it uses the exact same IAU 2000B, EOP and polar-motion transform as the
 * prediction engine. It is built as a standalone tool, not linked into the
 * flight service.
 */

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/calgnss_orbit_mini.c"

static void ecef_to_geodetic(
    cg_vec3_t ecef, double *latitude_deg, double *longitude_deg, double *height_m)
{
    const double a = 6378137.0;
    const double f = 1.0 / 298.257223563;
    const double e2 = f * (2.0 - f);
    const double p = hypot(ecef.x, ecef.y);
    double latitude = atan2(ecef.z, p * (1.0 - e2));
    double height = 0.0;
    int iteration;

    for (iteration = 0; iteration < 12; ++iteration) {
        const double sin_latitude = sin(latitude);
        const double n = a / sqrt(1.0 - e2 * sin_latitude * sin_latitude);
        const double next_height = p / cos(latitude) - n;
        const double next_latitude = atan2(
            ecef.z, p * (1.0 - e2 * n / (n + next_height)));
        height = next_height;
        if (fabs(next_latitude - latitude) < 1e-14) {
            latitude = next_latitude;
            break;
        }
        latitude = next_latitude;
    }
    *latitude_deg = latitude * 180.0 / CG_PI;
    *longitude_deg = atan2(ecef.y, ecef.x) * 180.0 / CG_PI;
    *height_m = height;
}

int main(int argc, char **argv)
{
    FILE *input;
    FILE *output;
    char line[4096];
    long long timestamp_ms;
    unsigned long frame_index;
    unsigned long point_index;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    size_t converted = 0U;

    if (argc != 3) {
        fprintf(stderr, "usage: %s orbit_points.csv locations.csv\n", argv[0]);
        return 2;
    }
    input = fopen(argv[1], "r");
    if (!input) {
        fprintf(stderr, "cannot open %s: %s\n", argv[1], strerror(errno));
        return 2;
    }
    output = fopen(argv[2], "w");
    if (!output) {
        fprintf(stderr, "cannot create %s: %s\n", argv[2], strerror(errno));
        fclose(input);
        return 2;
    }
    fprintf(output,
            "timestamp_ms,ecef_x_m,ecef_y_m,ecef_z_m,latitude_deg,longitude_deg,height_m\n");
    if (!fgets(line, sizeof(line), input)) {
        fprintf(stderr, "input CSV has no header\n");
        fclose(output);
        fclose(input);
        return 3;
    }
    while (fgets(line, sizeof(line), input)) {
        cg_time_t time_utc;
        cg_vec3_t j2000;
        cg_vec3_t ecef;
        double matrix[3][3];
        double latitude_deg;
        double longitude_deg;
        double height_m;
        int parsed = sscanf(
            line, "%lu,%lu,%lld,%lf,%lf,%lf,%lf,%lf,%lf",
            &frame_index, &point_index, &timestamp_ms,
            &x, &y, &z, &vx, &vy, &vz);
        (void)frame_index;
        (void)point_index;
        (void)vx;
        (void)vy;
        (void)vz;
        if (parsed != 9) {
            fprintf(stderr, "invalid input row: %s", line);
            fclose(output);
            fclose(input);
            return 3;
        }
        time_utc = cg_time_from_unix_seconds((double)timestamp_ms / 1000.0);
        cg_itrs_to_j2000_matrix(&time_utc, matrix);
        j2000.x = x;
        j2000.y = y;
        j2000.z = z;
        ecef.x = matrix[0][0] * j2000.x + matrix[1][0] * j2000.y +
                 matrix[2][0] * j2000.z;
        ecef.y = matrix[0][1] * j2000.x + matrix[1][1] * j2000.y +
                 matrix[2][1] * j2000.z;
        ecef.z = matrix[0][2] * j2000.x + matrix[1][2] * j2000.y +
                 matrix[2][2] * j2000.z;
        ecef_to_geodetic(ecef, &latitude_deg, &longitude_deg, &height_m);
        fprintf(output, "%lld,%.6f,%.6f,%.6f,%.10f,%.10f,%.4f\n",
                timestamp_ms, ecef.x, ecef.y, ecef.z,
                latitude_deg, longitude_deg, height_m);
        ++converted;
    }
    fclose(output);
    fclose(input);
    printf("converted=%zu output=%s\n", converted, argv[2]);
    return 0;
}
