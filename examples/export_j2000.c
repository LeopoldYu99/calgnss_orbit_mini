#include "calgnss.h"

#include <stdio.h>

static const char *k_input_csv = "..\\20260518\\67113_LLA_Position_2hour_step600s.csv";
static const char *k_output_csv = "..\\20260518\\67113_J2000_Calculated.csv";
static const char *k_start_time_text = "18 May 2026 04:00:00.000";
static const char *k_end_time_text = "18 May 2026 06:59:59.000";
static const double k_step_seconds = 1.0;

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

    status = cg_export_j2000_csv(input_csv, output_csv, start_ptr, end_ptr, step_seconds, &options);
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
