#include "calgnss_orbit_mini.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CSV_LINE_BUFFER 2048
#define CSV_MAX_FIELDS 32
#define EXAMPLE_STREAM_OBSERVATION_CAPACITY 600

static const char *k_input_csv = "E:\\JC\\GNSS_OrbitDetermination\\20260522\\67113_LLA.csv";
static const char *k_output_csv = "E:\\JC\\GNSS_OrbitDetermination\\20260522\\67113_J2000_Calculated.csv";
static const char *k_start_time_text = "18 May 2026 04:00:00.000";
static const char *k_end_time_text = "18 May 2026 06:00:00.000";
static const double k_step_seconds = 0.25;
static const double k_fit_window_minutes = 10;
static const int k_fit_degree = 9;

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

static int64_t example_days_from_civil(int y, unsigned m, unsigned d)
{
    int era;
    unsigned yoe;
    unsigned doy;
    unsigned doe;
    int mp;
    y -= (m <= 2U);
    era = (y >= 0 ? y : y - 399) / 400;
    yoe = (unsigned)(y - era * 400);
    mp = (int)m + ((m > 2U) ? -3 : 9);
    doy = (unsigned)((153 * mp + 2) / 5) + d - 1U;
    doe = yoe * 365U + yoe / 4U - yoe / 100U + doy;
    return (int64_t)era * 146097 + (int64_t)doe - 719468;
}

static void example_civil_from_days(int64_t z, int *y, unsigned *m, unsigned *d)
{
    int era;
    unsigned doe;
    unsigned yoe;
    unsigned doy;
    unsigned mp;
    z += 719468;
    era = (int)((z >= 0 ? z : z - 146096) / 146097);
    doe = (unsigned)(z - (int64_t)era * 146097);
    yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    *y = (int)yoe + era * 400;
    doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    mp = (5 * doy + 2) / 153;
    *d = doy - (153 * mp + 2) / 5 + 1;
    *m = mp + (mp < 10 ? 3U : (unsigned)-9);
    *y += (*m <= 2U);
}

static double example_jd_from_calendar(int y, int m, int d, int hour, int minute, double second)
{
    int yy = y;
    int mm = m;
    int a;
    int b;
    double day;
    if (mm <= 2) {
        yy -= 1;
        mm += 12;
    }
    a = yy / 100;
    b = 2 - a + a / 4;
    day = (double)d + ((double)hour + ((double)minute + second / 60.0) / 60.0) / 24.0;
    return floor(365.25 * (double)(yy + 4716)) +
           floor(30.6001 * (double)(mm + 1)) +
           day + (double)b - 1524.5;
}

static cg_time_t example_time_from_calendar(int y, int m, int d, int hour, int minute, double second)
{
    cg_time_t t;
    int whole_second = (int)floor(second);
    double frac = second - (double)whole_second;
    int64_t days = example_days_from_civil(y, (unsigned)m, (unsigned)d);
    t.year = y;
    t.month = m;
    t.day = d;
    t.hour = hour;
    t.minute = minute;
    t.second = second;
    t.jd_utc = example_jd_from_calendar(y, m, d, hour, minute, second);
    t.unix_seconds = (double)days * 86400.0 + (double)hour * 3600.0 +
                     (double)minute * 60.0 + (double)whole_second + frac;
    return t;
}

static cg_time_t example_time_add_seconds(const cg_time_t *time_utc, double seconds)
{
    double u = time_utc->unix_seconds + seconds;
    double whole_d = floor(u);
    int64_t whole = (int64_t)whole_d;
    double frac = u - whole_d;
    int64_t days = whole / 86400;
    int64_t sod = whole % 86400;
    int y;
    unsigned m;
    unsigned d;
    int hour;
    int minute;
    if (sod < 0) {
        sod += 86400;
        days -= 1;
    }
    example_civil_from_days(days, &y, &m, &d);
    hour = (int)(sod / 3600);
    minute = (int)((sod % 3600) / 60);
    return example_time_from_calendar(y, (int)m, (int)d, hour, minute,
                                      (double)(sod % 60) + frac);
}

static int example_month_number(const char *month)
{
    static const char *names[] = {
        "jan", "feb", "mar", "apr", "may", "jun",
        "jul", "aug", "sep", "oct", "nov", "dec"
    };
    char lower[4];
    int i;
    for (i = 0; i < 3 && month[i]; ++i) {
        lower[i] = (char)tolower((unsigned char)month[i]);
    }
    lower[i] = '\0';
    for (i = 0; i < 12; ++i) {
        if (strcmp(lower, names[i]) == 0) {
            return i + 1;
        }
    }
    return 0;
}

static int example_parse_time(const char *text, cg_time_t *out_time)
{
    char buf[128];
    char mon[16];
    char sep = 0;
    int y = 0, m = 0, d = 0, hh = 0, mm = 0;
    double ss = 0.0;
    size_t n;
    char *s;

    if (!text || !out_time) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    n = strlen(text);
    if (n >= sizeof(buf)) {
        return CG_ERR_PARSE;
    }
    memcpy(buf, text, n + 1);
    s = csv_trim(buf);
    n = strlen(s);
    if (n > 0 && s[n - 1] == 'Z') {
        s[n - 1] = '\0';
    }

    if (sscanf(s, "%d %15s %d %d:%d:%lf", &d, mon, &y, &hh, &mm, &ss) == 6) {
        m = example_month_number(mon);
        if (m > 0) {
            *out_time = example_time_from_calendar(y, m, d, hh, mm, ss);
            return CG_OK;
        }
    }

    if (sscanf(s, "%d-%d-%d%c%d:%d:%lf", &y, &m, &d, &sep, &hh, &mm, &ss) == 7 &&
        (sep == ' ' || sep == 'T')) {
        *out_time = example_time_from_calendar(y, m, d, hh, mm, ss);
        return CG_OK;
    }

    if (sscanf(s, "%d/%d/%d%c%d:%d:%lf", &y, &m, &d, &sep, &hh, &mm, &ss) == 7 &&
        (sep == ' ' || sep == 'T')) {
        *out_time = example_time_from_calendar(y, m, d, hh, mm, ss);
        return CG_OK;
    }

    return CG_ERR_PARSE;
}

static void example_format_time_iso(const cg_time_t *time_utc, char *buffer, size_t buffer_size)
{
    double unix_ms_d;
    int64_t unix_ms;
    int64_t whole_seconds;
    int milliseconds;
    int64_t days;
    int64_t sod;
    int y;
    unsigned m;
    unsigned d;
    int hour;
    int minute;
    int second;

    if (!time_utc || !buffer || buffer_size == 0) {
        return;
    }

    unix_ms_d = floor(time_utc->unix_seconds * 1000.0 + 0.5);
    unix_ms = (int64_t)unix_ms_d;
    whole_seconds = unix_ms / 1000;
    milliseconds = (int)(unix_ms % 1000);
    if (milliseconds < 0) {
        milliseconds += 1000;
        whole_seconds -= 1;
    }

    days = whole_seconds / 86400;
    sod = whole_seconds % 86400;
    if (sod < 0) {
        sod += 86400;
        days -= 1;
    }

    example_civil_from_days(days, &y, &m, &d);
    hour = (int)(sod / 3600);
    minute = (int)((sod % 3600) / 60);
    second = (int)(sod % 60);

    snprintf(buffer, buffer_size, "%04d-%02u-%02u %02d:%02d:%02d.%03d",
             y, m, d, hour, minute, second, milliseconds);
}

static int csv_detect_delimiter(const char *line)
{
    int comma_count = 0;
    int semicolon_count = 0;
    int in_quotes = 0;
    const char *p = line;

    while (*p) {
        if (*p == '"') {
            in_quotes = !in_quotes;
        } else if (!in_quotes) {
            if (*p == ',') {
                ++comma_count;
            } else if (*p == ';') {
                ++semicolon_count;
            } else if (*p == '\r' || *p == '\n') {
                break;
            }
        }
        ++p;
    }

    return semicolon_count > comma_count ? ';' : ',';
}

static int csv_split_line(char *line, char **fields, int max_fields, int delimiter)
{
    int count = 0;
    int in_quotes = 0;
    char *p = line;
    char *start = line;
    while (*p) {
        if (*p == '"') {
            in_quotes = !in_quotes;
        } else if (*p == delimiter && !in_quotes) {
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

typedef struct csv_lla_columns_t {
    int time_col;
    int lat_col;
    int lon_col;
    int alt_col;
    int alt_is_km;
    int delimiter;
} csv_lla_columns_t;

static int csv_read_lla_header(FILE *f, csv_lla_columns_t *columns)
{
    char line[CSV_LINE_BUFFER];
    char normalized_header[CSV_LINE_BUFFER];
    char *fields[CSV_MAX_FIELDS];
    int field_count;

    if (!f || !columns) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    if (!fgets(line, sizeof(line), f)) {
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

    columns->delimiter = csv_detect_delimiter(normalized_header);
    field_count = csv_split_line(normalized_header, fields, CSV_MAX_FIELDS, columns->delimiter);
    columns->time_col = csv_find_field(fields, field_count, "Time (UTCG)");
    columns->lat_col = csv_find_field(fields, field_count, "Lat (deg)");
    columns->lon_col = csv_find_field(fields, field_count, "Lon (deg)");
    columns->alt_col = csv_find_field(fields, field_count, "Alt (km)");
    if (columns->alt_col < 0) {
        columns->alt_col = csv_find_field_contains(fields, field_count, "Alt");
    }
    if (columns->time_col < 0 || columns->lat_col < 0 ||
        columns->lon_col < 0 || columns->alt_col < 0) {
        return CG_ERR_PARSE;
    }
    columns->alt_is_km =
        strstr(fields[columns->alt_col], "(km)") != NULL ||
        strstr(fields[columns->alt_col], "(KM)") != NULL;
    return CG_OK;
}

static int csv_parse_lla_observation(
    char *line,
    const csv_lla_columns_t *columns,
    cg_observation_t *out_observation,
    int *out_skip)
{
    char *fields[CSV_MAX_FIELDS];
    int field_count;
    int max_col;
    char *endptr;
    double alt_value;

    if (!line || !columns || !out_observation || !out_skip) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    *out_skip = 0;
    max_col = columns->time_col;
    if (columns->lat_col > max_col) max_col = columns->lat_col;
    if (columns->lon_col > max_col) max_col = columns->lon_col;
    if (columns->alt_col > max_col) max_col = columns->alt_col;

    field_count = csv_split_line(line, fields, CSV_MAX_FIELDS, columns->delimiter);
    if (field_count <= max_col) {
        *out_skip = 1;
        return CG_OK;
    }

    memset(out_observation, 0, sizeof(*out_observation));
    if (example_parse_time(fields[columns->time_col], &out_observation->time_utc) != CG_OK) {
        return CG_ERR_PARSE;
    }

    errno = 0;
    out_observation->lat_deg = strtod(fields[columns->lat_col], &endptr);
    if (errno || endptr == fields[columns->lat_col]) {
        return CG_ERR_PARSE;
    }

    errno = 0;
    out_observation->lon_deg = strtod(fields[columns->lon_col], &endptr);
    if (errno || endptr == fields[columns->lon_col]) {
        return CG_ERR_PARSE;
    }

    errno = 0;
    alt_value = strtod(fields[columns->alt_col], &endptr);
    if (errno || endptr == fields[columns->alt_col]) {
        return CG_ERR_PARSE;
    }
    out_observation->alt_m = columns->alt_is_km ? alt_value * 1000.0 : alt_value;
    return CG_OK;
}

static void write_j2000_csv_row(FILE *out, const cg_time_t *time_utc, const cg_state_t *state)
{
    char time_text[64];
    example_format_time_iso(time_utc, time_text, sizeof(time_text));
    fprintf(out, "%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            time_text,
            state->r_j2000_m.x / 1000.0,
            state->r_j2000_m.y / 1000.0,
            state->r_j2000_m.z / 1000.0,
            state->v_j2000_mps.x / 1000.0,
            state->v_j2000_mps.y / 1000.0,
            state->v_j2000_mps.z / 1000.0);
}

static int export_j2000_csv(
    const char *input_lla_csv,
    const char *output_j2000_csv,
    const cg_time_t *start_time_utc,
    const cg_time_t *end_time_utc,
    double step_seconds,
    const cg_options_t *options)
{
    cg_observation_t stream_buffer[EXAMPLE_STREAM_OBSERVATION_CAPACITY];
    cg_context_t *context = NULL;
    csv_lla_columns_t columns;
    char line[CSV_LINE_BUFFER];
    cg_time_t current;
    FILE *in;
    FILE *out;
    size_t input_count = 0;
    int have_current = 0;
    int have_observation = 0;
    int status;

    if (!input_lla_csv || !output_j2000_csv || step_seconds <= 0.0) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    memset(&current, 0, sizeof(current));

    in = fopen(input_lla_csv, "rb");
    if (!in) {
        return CG_ERR_IO;
    }

    status = csv_read_lla_header(in, &columns);
    if (status != CG_OK) {
        fclose(in);
        return status;
    }

    status = cg_context_create(&context, stream_buffer, EXAMPLE_STREAM_OBSERVATION_CAPACITY, options);
    if (status != CG_OK) {
        fclose(in);
        return status;
    }

    out = fopen(output_j2000_csv, "wb");
    if (!out) {
        fclose(in);
        cg_context_destroy(context);
        return CG_ERR_IO;
    }

    fprintf(out, "\"Time (UTCG)\",\"x (km)\",\"y (km)\",\"z (km)\",\"vx (km/sec)\",\"vy (km/sec)\",\"vz (km/sec)\"\n");

    while (fgets(line, sizeof(line), in)) {
        cg_observation_t observation;
        int skip = 0;

        status = csv_parse_lla_observation(line, &columns, &observation, &skip);
        if (status != CG_OK) {
            fclose(in);
            fclose(out);
            cg_context_destroy(context);
            return status;
        }
        if (skip) {
            continue;
        }

        if (!have_current) {
            current = start_time_utc ? *start_time_utc : observation.time_utc;
            have_current = 1;
        }

        status = cg_context_push(context, &observation);
        if (status != CG_OK) {
            fclose(in);
            fclose(out);
            cg_context_destroy(context);
            return status;
        }
        ++input_count;
        have_observation = 1;

        while (cg_context_count(context) >= 2 &&
               current.unix_seconds <= observation.time_utc.unix_seconds + 1.0e-9 &&
               (!end_time_utc || current.unix_seconds <= end_time_utc->unix_seconds + 1.0e-9)) {
            cg_state_t state;

            status = cg_context_query_state(context, &current, &state);
            if (status != CG_OK) {
                fclose(in);
                fclose(out);
                cg_context_destroy(context);
                return status;
            }

            write_j2000_csv_row(out, &current, &state);
            current = example_time_add_seconds(&current, step_seconds);
        }

        if (end_time_utc && current.unix_seconds > end_time_utc->unix_seconds + 1.0e-9) {
            break;
        }
    }

    if (!have_observation || input_count < 2) {
        fclose(in);
        fclose(out);
        cg_context_destroy(context);
        return CG_ERR_PARSE;
    }

    while (have_current &&
           end_time_utc &&
           cg_context_count(context) >= 2 &&
           current.unix_seconds <= end_time_utc->unix_seconds + 1.0e-9) {
        cg_state_t state;

        status = cg_context_query_state(context, &current, &state);
        if (status != CG_OK) {
            fclose(in);
            fclose(out);
            cg_context_destroy(context);
            return status;
        }

        write_j2000_csv_row(out, &current, &state);
        current = example_time_add_seconds(&current, step_seconds);
    }

    fclose(in);
    fclose(out);
    cg_context_destroy(context);
    return CG_OK;
}

static int parse_config_time(const char *label, const char *text, cg_time_t *time, cg_time_t **time_ptr)
{
    int status;
    if (!text || text[0] == '\0') {
        *time_ptr = NULL;
        return CG_OK;
    }

    status = example_parse_time(text, time);
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

    options.fit_window_minutes = k_fit_window_minutes;
    options.degree = k_fit_degree;

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
    fprintf(stderr, "Observation capacity: %d\n", EXAMPLE_STREAM_OBSERVATION_CAPACITY);
    fprintf(stderr, "Fit window: %.3f min, degree %d\n", options.fit_window_minutes, options.degree);
    fprintf(stderr, "Extrapolation history limit: %.3f s\n", options.extrapolation_history_seconds);
    fprintf(stderr, "Future model: J2 RK4 + previous-orbit residual correction\n");
    return 0;
}
