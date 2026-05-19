#include "calgnss.h"

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CSV_LINE_BUFFER 2048
#define CSV_MAX_FIELDS 32

static const char *k_input_csv = "..\\20260518\\67217_LLA_Position_2hour.csv";
static const char *k_output_csv = "..\\20260518\\67217_J2000_Calculated.csv";
static const char *k_start_time_text = "18 May 2026 04:00:00.000";
static const char *k_end_time_text = "18 May 2026 06:59:59.000";
static const double k_step_seconds = 1.0;

static char *csv_trim(char *s)
{
    char *end;
    while (*s && isspace((unsigned char)*s)) {
        ++s;
    }
    end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) {
        *--end = '\0';
    }
    if (*s == '"' && end > s + 1 && end[-1] == '"') {
        *--end = '\0';
        ++s;
    }
    return s;
}

static int csv_split_line(char *line, char **fields, int max_fields)
{
    int count = 0;
    int in_quotes = 0;
    char *p = line;
    char *start = line;
    while (*p) {
        if (*p == '"') {
            in_quotes = !in_quotes;
        } else if (*p == ',' && !in_quotes) {
            *p = '\0';
            if (count < max_fields) {
                fields[count++] = csv_trim(start);
            }
            start = p + 1;
        } else if (*p == '\r' || *p == '\n') {
            *p = '\0';
            break;
        }
        ++p;
    }
    if (count < max_fields) {
        fields[count++] = csv_trim(start);
    }
    return count;
}

static int csv_field_equals(const char *a, const char *b)
{
    while (*a && *b) {
        if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) {
            return 0;
        }
        ++a;
        ++b;
    }
    return *a == '\0' && *b == '\0';
}

static int csv_find_field(char **fields, int count, const char *name)
{
    int i;
    for (i = 0; i < count; ++i) {
        if (csv_field_equals(fields[i], name)) {
            return i;
        }
    }
    return -1;
}

static int csv_find_field_contains(char **fields, int count, const char *needle)
{
    int i;
    for (i = 0; i < count; ++i) {
        if (strstr(fields[i], needle) != NULL) {
            return i;
        }
    }
    return -1;
}

static int compare_observations(const void *a, const void *b)
{
    const cg_observation_t *oa = (const cg_observation_t *)a;
    const cg_observation_t *ob = (const cg_observation_t *)b;
    if (oa->time_utc.unix_seconds < ob->time_utc.unix_seconds) {
        return -1;
    }
    if (oa->time_utc.unix_seconds > ob->time_utc.unix_seconds) {
        return 1;
    }
    return 0;
}

static int load_lla_csv(const char *path, cg_observation_t **out_observations, size_t *out_count)
{
    FILE *f;
    char line[CSV_LINE_BUFFER];
    char normalized_header[CSV_LINE_BUFFER];
    char *fields[CSV_MAX_FIELDS];
    int field_count;
    int time_col;
    int lat_col;
    int lon_col;
    int alt_col;
    int lat_rate_col;
    int lon_rate_col;
    int alt_rate_col;
    int alt_is_km;
    cg_observation_t *observations = NULL;
    size_t count = 0;
    size_t capacity = 0;

    if (!path || !out_observations || !out_count) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    *out_observations = NULL;
    *out_count = 0;

    f = fopen(path, "rb");
    if (!f) {
        return CG_ERR_IO;
    }
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        return CG_ERR_PARSE;
    }
    {
        size_t si = 0;
        size_t di = 0;
        if ((unsigned char)line[0] == 0xEF &&
            (unsigned char)line[1] == 0xBB &&
            (unsigned char)line[2] == 0xBF) {
            si = 3;
        }
        while (line[si] && di + 1 < sizeof(normalized_header)) {
            if (line[si] != '"') {
                normalized_header[di++] = line[si];
            }
            ++si;
        }
        normalized_header[di] = '\0';
    }
    field_count = csv_split_line(normalized_header, fields, CSV_MAX_FIELDS);
    time_col = csv_find_field(fields, field_count, "Time (UTCG)");
    lat_col = csv_find_field(fields, field_count, "Lat (deg)");
    lon_col = csv_find_field(fields, field_count, "Lon (deg)");
    alt_col = csv_find_field(fields, field_count, "Alt (km)");
    if (alt_col < 0) {
        alt_col = csv_find_field_contains(fields, field_count, "Alt");
    }
    lat_rate_col = csv_find_field(fields, field_count, "Lat Rate (deg/sec)");
    lon_rate_col = csv_find_field(fields, field_count, "Lon Rate (deg/sec)");
    alt_rate_col = csv_find_field(fields, field_count, "Alt Rate (km/sec)");
    if (time_col < 0 || lat_col < 0 || lon_col < 0 || alt_col < 0) {
        fclose(f);
        return CG_ERR_PARSE;
    }
    alt_is_km = strstr(fields[alt_col], "(km)") != NULL || strstr(fields[alt_col], "(KM)") != NULL;

    while (fgets(line, sizeof(line), f)) {
        cg_observation_t item;
        char *endptr;
        double alt_value;
        field_count = csv_split_line(line, fields, CSV_MAX_FIELDS);
        if (field_count <= alt_col) {
            continue;
        }
        memset(&item, 0, sizeof(item));
        if (cg_parse_time(fields[time_col], &item.time_utc) != CG_OK) {
            free(observations);
            fclose(f);
            return CG_ERR_PARSE;
        }
        errno = 0;
        item.lat_deg = strtod(fields[lat_col], &endptr);
        if (errno || endptr == fields[lat_col]) {
            free(observations);
            fclose(f);
            return CG_ERR_PARSE;
        }
        item.lon_deg = strtod(fields[lon_col], &endptr);
        alt_value = strtod(fields[alt_col], &endptr);
        item.alt_m = alt_is_km ? alt_value * 1000.0 : alt_value;
        item.has_rates = lat_rate_col >= 0 && lon_rate_col >= 0 && alt_rate_col >= 0 &&
                         field_count > lat_rate_col && field_count > lon_rate_col && field_count > alt_rate_col;
        if (item.has_rates) {
            item.lat_rate_degps = strtod(fields[lat_rate_col], NULL);
            item.lon_rate_degps = strtod(fields[lon_rate_col], NULL);
            item.alt_rate_mps = strtod(fields[alt_rate_col], NULL) * 1000.0;
        }
        if (count == capacity) {
            size_t new_capacity = capacity == 0 ? 128 : capacity * 2;
            cg_observation_t *new_observations =
                (cg_observation_t *)realloc(observations, new_capacity * sizeof(*observations));
            if (!new_observations) {
                free(observations);
                fclose(f);
                return CG_ERR_NO_MEMORY;
            }
            observations = new_observations;
            capacity = new_capacity;
        }
        observations[count++] = item;
    }
    fclose(f);
    if (count < 2) {
        free(observations);
        return CG_ERR_PARSE;
    }

    qsort(observations, count, sizeof(*observations), compare_observations);
    *out_observations = observations;
    *out_count = count;
    return CG_OK;
}

static int export_j2000_csv(
    const char *input_lla_csv,
    const char *output_j2000_csv,
    const cg_time_t *start_time_utc,
    const cg_time_t *end_time_utc,
    double step_seconds,
    const cg_options_t *options)
{
    cg_observation_t *observations = NULL;
    size_t count = 0;
    cg_time_t start;
    cg_time_t end;
    cg_time_t current;
    FILE *out;
    int status;

    if (!input_lla_csv || !output_j2000_csv || step_seconds <= 0.0) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    status = load_lla_csv(input_lla_csv, &observations, &count);
    if (status != CG_OK) {
        return status;
    }

    status = cg_precompute_observations(observations, count);
    if (status != CG_OK) {
        free(observations);
        return status;
    }

    start = start_time_utc ? *start_time_utc : observations[0].time_utc;
    end = end_time_utc ? *end_time_utc : observations[count - 1].time_utc;

    out = fopen(output_j2000_csv, "wb");
    if (!out) {
        free(observations);
        return CG_ERR_IO;
    }

    fprintf(out, "\"Time (UTCG)\",\"x (km)\",\"y (km)\",\"z (km)\",\"vx (km/sec)\",\"vy (km/sec)\",\"vz (km/sec)\"\n");
    current = start;
    while (current.unix_seconds <= end.unix_seconds + 1.0e-9) {
        cg_state_t state;
        char time_text[64];

        status = cg_query_state(observations, count, &current, options, &state);
        if (status != CG_OK) {
            fclose(out);
            free(observations);
            return status;
        }

        cg_format_time_iso(&current, time_text, sizeof(time_text));
        fprintf(out, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                time_text,
                state.r_j2000_m.x / 1000.0,
                state.r_j2000_m.y / 1000.0,
                state.r_j2000_m.z / 1000.0,
                state.v_j2000_mps.x / 1000.0,
                state.v_j2000_mps.y / 1000.0,
                state.v_j2000_mps.z / 1000.0);
        current = cg_time_add_seconds(&current, step_seconds);
    }

    fclose(out);
    free(observations);
    return CG_OK;
}

static int parse_config_time(const char *label, const char *text, cg_time_t *time, cg_time_t **time_ptr)
{
    int status;
    if (!text || text[0] == '\0') {
        *time_ptr = NULL;
        return CG_OK;
    }

    status = cg_parse_time(text, time);
    if (status != CG_OK) {
        fprintf(stderr, "Invalid %s time: %s\n", label, text);
        return status;
    }

    *time_ptr = time;
    return CG_OK;
}

int main(void)
{
    const char *input_csv = k_input_csv;
    const char *output_csv = k_output_csv;
    cg_time_t start_time;
    cg_time_t end_time;
    cg_time_t *start_ptr = NULL;
    cg_time_t *end_ptr = NULL;
    cg_options_t options = cg_default_options();
    double step_seconds = k_step_seconds;
    int status;

    status = parse_config_time("start", k_start_time_text, &start_time, &start_ptr);
    if (status != CG_OK) {
        return 2;
    }
    status = parse_config_time("end", k_end_time_text, &end_time, &end_ptr);
    if (status != CG_OK) {
        return 2;
    }

    status = export_j2000_csv(input_csv, output_csv, start_ptr, end_ptr, step_seconds, &options);
    if (status != CG_OK) {
        fprintf(stderr, "Export failed: %s (%d)\n", cg_status_string(status), status);
        fprintf(stderr, "Input:  %s\nOutput: %s\n", input_csv, output_csv);
        return 1;
    }
    fprintf(stderr, "Input CSV:  %s\n", input_csv);
    fprintf(stderr, "Output CSV: %s\n", output_csv);
    fprintf(stderr, "Step: %.3f s\n", step_seconds);
    fprintf(stderr, "Fit window: %.3f min, degree %d\n", options.fit_window_minutes, options.degree);
    fprintf(stderr, "Future model: J2 RK4 + previous-orbit residual correction\n");
    return 0;
}
