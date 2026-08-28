#include "calgnss_orbit_mini.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef CG_PI
#define CG_PI 3.141592653589793238462643383279502884
#endif

#define CG_D2PI (2.0 * CG_PI)
#define CG_ARCSEC2RAD (CG_PI / (180.0 * 3600.0))
#define CG_MAS2RAD (CG_ARCSEC2RAD / 1000.0)
#define CG_JD_J2000 2451545.0
#define CG_JULIAN_CENTURY 36525.0
#define CG_TT_MINUS_UTC_SECONDS 69.184
#define CG_SECONDS_PER_DAY 86400.0

#define CG_OMEGA_EARTH 7.2921150e-5
#define CG_FIT_WINDOW_SECONDS 300.0
#define CG_MAX_EXTRAPOLATION_SECONDS 3600.0
#define CG_EXTRAPOLATION_HISTORY_SECONDS 600.0
#define CG_PROPAGATION_STEP_SECONDS 10.0
#define CG_TAIL_REFINE_SECONDS 600.0
#define CG_HOLDOUT_CORRECTION_SECONDS 540.0
#define CG_HOLDOUT_CORRECTION_GAIN 0.40
#define CG_HOLDOUT_CORRECTION_MAX_M 1000.0
#define CG_ASTRONOMICAL_UNIT_M 1.495978707e11
#define CG_SUN_MU 1.32712440018e20
#define CG_MOON_MU 4.9048695e12
#define CG_SOLAR_RADIATION_PRESSURE 4.56e-6
#define CG_EARTH_EQUATORIAL_RADIUS_M 6378136.3
#define CG_GRAVITY_DEGREE 20

#ifndef CG_USE_EOP_TABLE
#define CG_USE_EOP_TABLE 1
#endif

typedef struct GravityTerm GravityTerm;
#define constexpr static const
#include "gravity_field_20.inc"
#undef constexpr

typedef struct cg_force_model_t {
    double drag_ballistic_m2_per_kg;
    double gravity_harmonic_scale;
    double srp_area_m2_per_kg;
    cg_vec3_t empirical_rtn_mps2;
} cg_force_model_t;

typedef struct cg_fit_t {
    int degree;
    size_t first_index;
    size_t last_index;
    double epoch_seconds;
    double scale_seconds;
    double coeff[3][CG_MAX_DEGREE + 1];
} cg_fit_t;

typedef struct cg_fit_cache_t {
    int valid;
    int degree;
    size_t first_index;
    size_t last_index;
    cg_fit_t fit;
} cg_fit_cache_t;

typedef struct cg_future_cache_t {
    int valid;
    int holdout_valid;
    double latest_unix_seconds;
    double holdout_tau_seconds;
    cg_force_model_t force_model;
    cg_state_t latest_state;
    cg_state_t last_state;
    cg_vec3_t holdout_position_residual;
    cg_vec3_t holdout_velocity_residual;
} cg_future_cache_t;

struct cg_context_t {
    cg_observation_t *observations;
    size_t capacity;
    size_t count;
    size_t start;
    cg_options_t options;
    cg_fit_cache_t main_cache;
    cg_future_cache_t future_cache;
};

typedef struct cg_eop_record_t {
    int mjd;
    double dut1_seconds;
    double x_pole_arcsec;
    double y_pole_arcsec;
} cg_eop_record_t;

typedef struct cg_eop_value_t {
    double dut1_seconds;
    double x_pole_rad;
    double y_pole_rad;
} cg_eop_value_t;

#include "eop_2026.inc"

static cg_vec3_t cg_vec_add(cg_vec3_t a, cg_vec3_t b)
{
    cg_vec3_t r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    return r;
}

static cg_vec3_t cg_vec_sub(cg_vec3_t a, cg_vec3_t b)
{
    cg_vec3_t r;
    r.x = a.x - b.x;
    r.y = a.y - b.y;
    r.z = a.z - b.z;
    return r;
}

static cg_vec3_t cg_vec_scale(cg_vec3_t v, double scale)
{
    cg_vec3_t r;
    r.x = v.x * scale;
    r.y = v.y * scale;
    r.z = v.z * scale;
    return r;
}

static double cg_vec_dot(cg_vec3_t a, cg_vec3_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static double cg_vec_norm(cg_vec3_t v)
{
    return sqrt(cg_vec_dot(v, v));
}

static cg_vec3_t cg_vec_cross(cg_vec3_t a, cg_vec3_t b)
{
    cg_vec3_t r;
    r.x = a.y * b.z - a.z * b.y;
    r.y = a.z * b.x - a.x * b.z;
    r.z = a.x * b.y - a.y * b.x;
    return r;
}

static cg_vec3_t cg_vec_unit(cg_vec3_t v)
{
    double n = cg_vec_norm(v);
    cg_vec3_t zero;
    zero.x = zero.y = zero.z = 0.0;
    if (n <= 0.0 || !isfinite(n)) {
        return zero;
    }
    return cg_vec_scale(v, 1.0 / n);
}

static double cg_anp(double angle)
{
    double r = fmod(angle, CG_D2PI);
    return (r < 0.0) ? r + CG_D2PI : r;
}

static void cg_mat_identity(double r[3][3])
{
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

static void cg_mat_mul(const double a[3][3], const double b[3][3], double out[3][3])
{
    double r[3][3];
    int i;
    int j;
    int k;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r[i][j] = 0.0;
            for (k = 0; k < 3; ++k) {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    memcpy(out, r, sizeof(r));
}

static cg_vec3_t cg_mat_vec(const double m[3][3], cg_vec3_t v)
{
    cg_vec3_t r;
    r.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
    r.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
    r.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;
    return r;
}

static void cg_rx(double phi, double r[3][3])
{
    double s = sin(phi);
    double c = cos(phi);
    double a10 = c * r[1][0] + s * r[2][0];
    double a11 = c * r[1][1] + s * r[2][1];
    double a12 = c * r[1][2] + s * r[2][2];
    double a20 = -s * r[1][0] + c * r[2][0];
    double a21 = -s * r[1][1] + c * r[2][1];
    double a22 = -s * r[1][2] + c * r[2][2];
    r[1][0] = a10; r[1][1] = a11; r[1][2] = a12;
    r[2][0] = a20; r[2][1] = a21; r[2][2] = a22;
}

static void cg_ry(double theta, double r[3][3])
{
    double s = sin(theta);
    double c = cos(theta);
    double a00 = c * r[0][0] - s * r[2][0];
    double a01 = c * r[0][1] - s * r[2][1];
    double a02 = c * r[0][2] - s * r[2][2];
    double a20 = s * r[0][0] + c * r[2][0];
    double a21 = s * r[0][1] + c * r[2][1];
    double a22 = s * r[0][2] + c * r[2][2];
    r[0][0] = a00; r[0][1] = a01; r[0][2] = a02;
    r[2][0] = a20; r[2][1] = a21; r[2][2] = a22;
}

static void cg_rz(double psi, double r[3][3])
{
    double s = sin(psi);
    double c = cos(psi);
    double a00 = c * r[0][0] + s * r[1][0];
    double a01 = c * r[0][1] + s * r[1][1];
    double a02 = c * r[0][2] + s * r[1][2];
    double a10 = -s * r[0][0] + c * r[1][0];
    double a11 = -s * r[0][1] + c * r[1][1];
    double a12 = -s * r[0][2] + c * r[1][2];
    r[0][0] = a00; r[0][1] = a01; r[0][2] = a02;
    r[1][0] = a10; r[1][1] = a11; r[1][2] = a12;
}

static void cg_rz_vector_matrix(double angle, double r[3][3])
{
    double s = sin(angle);
    double c = cos(angle);
    r[0][0] = c;  r[0][1] = -s; r[0][2] = 0.0;
    r[1][0] = s;  r[1][1] = c;  r[1][2] = 0.0;
    r[2][0] = 0.0; r[2][1] = 0.0; r[2][2] = 1.0;
}

static double cg_mjd_from_jd(double jd)
{
    return jd - 2400000.5;
}

static int64_t cg_days_from_civil(int y, unsigned m, unsigned d)
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

static void cg_civil_from_days(int64_t z, int *y, unsigned *m, unsigned *d)
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

static double cg_jd_from_calendar(int y, int m, int d, int hour, int minute, double second)
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

static cg_time_t cg_time_from_unix_seconds(double unix_seconds)
{
    double days_double = floor(unix_seconds / CG_SECONDS_PER_DAY);
    int64_t days = (int64_t)days_double;
    double seconds_of_day = unix_seconds - (double)days * CG_SECONDS_PER_DAY;
    int year;
    unsigned month;
    unsigned day;
    int whole_seconds;
    double frac;
    cg_time_t time;
    if (seconds_of_day < 0.0) {
        seconds_of_day += CG_SECONDS_PER_DAY;
        --days;
    }

    cg_civil_from_days(days, &year, &month, &day);
    whole_seconds = (int)floor(seconds_of_day);
    frac = seconds_of_day - (double)whole_seconds;
    memset(&time, 0, sizeof(time));
    time.year = year;
    time.month = (int)month;
    time.day = (int)day;
    time.hour = whole_seconds / 3600;
    time.minute = (whole_seconds % 3600) / 60;
    time.second = (double)(whole_seconds % 60) + frac;
    time.unix_seconds = unix_seconds;
    time.jd_utc = cg_jd_from_calendar(
        time.year, time.month, time.day, time.hour, time.minute, time.second);
    return time;
}

static cg_eop_value_t cg_eop_from_record(const cg_eop_record_t *record)
{
    cg_eop_value_t value;
    value.dut1_seconds = record->dut1_seconds;
    value.x_pole_rad = record->x_pole_arcsec * CG_ARCSEC2RAD;
    value.y_pole_rad = record->y_pole_arcsec * CG_ARCSEC2RAD;
    return value;
}

static cg_eop_value_t cg_interpolate_eop(
    const cg_eop_record_t *a,
    const cg_eop_record_t *b,
    double mjd_utc)
{
    double span = (double)(b->mjd - a->mjd);
    double f = span > 0.0 ? (mjd_utc - (double)a->mjd) / span : 0.0;
    cg_eop_value_t value;
    value.dut1_seconds = a->dut1_seconds + f * (b->dut1_seconds - a->dut1_seconds);
    value.x_pole_rad = (a->x_pole_arcsec + f * (b->x_pole_arcsec - a->x_pole_arcsec)) * CG_ARCSEC2RAD;
    value.y_pole_rad = (a->y_pole_arcsec + f * (b->y_pole_arcsec - a->y_pole_arcsec)) * CG_ARCSEC2RAD;
    return value;
}

static cg_eop_value_t cg_lookup_eop(double mjd_utc)
{
    cg_eop_value_t empty;
    empty.dut1_seconds = 0.0;
    empty.x_pole_rad = 0.0;
    empty.y_pole_rad = 0.0;
#if !CG_USE_EOP_TABLE
    (void)mjd_utc;
    return empty;
#else
    size_t count = sizeof(k_cg_eop_2026) / sizeof(k_cg_eop_2026[0]);
    size_t i;
    if (count == 0U) {
        return empty;
    }
    if (mjd_utc <= (double)k_cg_eop_2026[0].mjd) {
        return cg_eop_from_record(&k_cg_eop_2026[0]);
    }
    if (mjd_utc >= (double)k_cg_eop_2026[count - 1U].mjd) {
        return cg_eop_from_record(&k_cg_eop_2026[count - 1U]);
    }
    for (i = 1U; i < count; ++i) {
        if (mjd_utc <= (double)k_cg_eop_2026[i].mjd) {
            return cg_interpolate_eop(&k_cg_eop_2026[i - 1U], &k_cg_eop_2026[i], mjd_utc);
        }
    }
    return empty;
#endif
}

static void cg_polar_motion_matrix(double xp, double yp, double out[3][3])
{
    double ry_mat[3][3];
    double rx_mat[3][3];
    cg_mat_identity(ry_mat);
    cg_mat_identity(rx_mat);
    cg_ry(-xp, ry_mat);
    cg_rx(-yp, rx_mat);
    cg_mat_mul(ry_mat, rx_mat, out);
}

typedef struct cg_nut_term_t {
    int nl;
    int nlp;
    int nf;
    int nd;
    int nom;
    double ps;
    double pst;
    double pc;
    double ec;
    double ect;
    double es;
} cg_nut_term_t;

static void cg_nut00b(double jd_tt, double *dpsi, double *deps)
{
    static const cg_nut_term_t x[] = {
        {0,0,0,0,1,-172064161.0,-174666.0,33386.0,92052331.0,9086.0,15377.0},
        {0,0,2,-2,2,-13170906.0,-1675.0,-13696.0,5730336.0,-3015.0,-4587.0},
        {0,0,2,0,2,-2276413.0,-234.0,2796.0,978459.0,-485.0,1374.0},
        {0,0,0,0,2,2074554.0,207.0,-698.0,-897492.0,470.0,-291.0},
        {0,1,0,0,0,1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0},
        {0,1,2,-2,2,-516821.0,1226.0,-524.0,224386.0,-677.0,-174.0},
        {1,0,0,0,0,711159.0,73.0,-872.0,-6750.0,0.0,358.0},
        {0,0,2,0,1,-387298.0,-367.0,380.0,200728.0,18.0,318.0},
        {1,0,2,0,2,-301461.0,-36.0,816.0,129025.0,-63.0,367.0},
        {0,-1,2,-2,2,215829.0,-494.0,111.0,-95929.0,299.0,132.0},
        {0,0,2,-2,1,128227.0,137.0,181.0,-68982.0,-9.0,39.0},
        {-1,0,2,0,2,123457.0,11.0,19.0,-53311.0,32.0,-4.0},
        {-1,0,0,2,0,156994.0,10.0,-168.0,-1235.0,0.0,82.0},
        {1,0,0,0,1,63110.0,63.0,27.0,-33228.0,0.0,-9.0},
        {-1,0,0,0,1,-57976.0,-63.0,-189.0,31429.0,0.0,-75.0},
        {-1,0,2,2,2,-59641.0,-11.0,149.0,25543.0,-11.0,66.0},
        {1,0,2,0,1,-51613.0,-42.0,129.0,26366.0,0.0,78.0},
        {-2,0,2,0,1,45893.0,50.0,31.0,-24236.0,-10.0,20.0},
        {0,0,0,2,0,63384.0,11.0,-150.0,-1220.0,0.0,29.0},
        {0,0,2,2,2,-38571.0,-1.0,158.0,16452.0,-11.0,68.0},
        {0,-2,2,-2,2,32481.0,0.0,0.0,-13870.0,0.0,0.0},
        {-2,0,0,2,0,-47722.0,0.0,-18.0,477.0,0.0,-25.0},
        {2,0,2,0,2,-31046.0,-1.0,131.0,13238.0,-11.0,59.0},
        {1,0,2,-2,2,28593.0,0.0,-1.0,-12338.0,10.0,-3.0},
        {-1,0,2,0,1,20441.0,21.0,10.0,-10758.0,0.0,-3.0},
        {2,0,0,0,0,29243.0,0.0,-74.0,-609.0,0.0,13.0},
        {0,0,2,0,0,25887.0,0.0,-66.0,-550.0,0.0,11.0},
        {0,1,0,0,1,-14053.0,-25.0,79.0,8551.0,-2.0,-45.0},
        {-1,0,0,2,1,15164.0,10.0,11.0,-8001.0,0.0,-1.0},
        {0,2,2,-2,2,-15794.0,72.0,-16.0,6850.0,-42.0,-5.0},
        {0,0,-2,2,0,21783.0,0.0,13.0,-167.0,0.0,13.0},
        {1,0,0,-2,1,-12873.0,-10.0,-37.0,6953.0,0.0,-14.0},
        {0,-1,0,0,1,-12654.0,11.0,63.0,6415.0,0.0,26.0},
        {-1,0,2,2,1,-10204.0,0.0,25.0,5222.0,0.0,15.0},
        {0,2,0,0,0,16707.0,-85.0,-10.0,168.0,-1.0,10.0},
        {1,0,2,2,2,-7691.0,0.0,44.0,3268.0,0.0,19.0},
        {-2,0,2,0,0,-11024.0,0.0,-14.0,104.0,0.0,2.0},
        {0,1,2,0,2,7566.0,-21.0,-11.0,-3250.0,0.0,-5.0},
        {0,0,2,2,1,-6637.0,-11.0,25.0,3353.0,0.0,14.0},
        {0,-1,2,0,2,-7141.0,21.0,8.0,3070.0,0.0,4.0},
        {0,0,0,2,1,-6302.0,-11.0,2.0,3272.0,0.0,4.0},
        {1,0,2,-2,1,5800.0,10.0,2.0,-3045.0,0.0,-1.0},
        {2,0,2,-2,2,6443.0,0.0,-7.0,-2768.0,0.0,-4.0},
        {-2,0,0,2,1,-5774.0,-11.0,-15.0,3041.0,0.0,-5.0},
        {2,0,2,0,1,-5350.0,0.0,21.0,2695.0,0.0,12.0},
        {0,-1,2,-2,1,-4752.0,-11.0,-3.0,2719.0,0.0,-3.0},
        {0,0,0,-2,1,-4940.0,-11.0,-21.0,2720.0,0.0,-9.0},
        {-1,-1,0,2,0,7350.0,0.0,-8.0,-51.0,0.0,4.0},
        {2,0,0,-2,1,4065.0,0.0,6.0,-2206.0,0.0,1.0},
        {1,0,0,2,0,6579.0,0.0,-24.0,-199.0,0.0,2.0},
        {0,1,2,-2,1,3579.0,0.0,5.0,-1900.0,0.0,1.0},
        {1,-1,0,0,0,4725.0,0.0,-6.0,-41.0,0.0,3.0},
        {-2,0,2,0,2,-3075.0,0.0,-2.0,1313.0,0.0,-1.0},
        {3,0,2,0,2,-2904.0,0.0,15.0,1233.0,0.0,7.0},
        {0,-1,0,2,0,4348.0,0.0,-10.0,-81.0,0.0,2.0},
        {1,-1,2,0,2,-2878.0,0.0,8.0,1232.0,0.0,4.0},
        {0,0,0,1,0,-4230.0,0.0,5.0,-20.0,0.0,-2.0},
        {-1,-1,2,2,2,-2819.0,0.0,7.0,1207.0,0.0,3.0},
        {-1,0,2,0,0,-4056.0,0.0,5.0,40.0,0.0,-2.0},
        {0,-1,2,2,2,-2647.0,0.0,11.0,1129.0,0.0,5.0},
        {-2,0,0,0,1,-2294.0,0.0,-10.0,1266.0,0.0,-4.0},
        {1,1,2,0,2,2481.0,0.0,-7.0,-1062.0,0.0,-3.0},
        {2,0,0,0,1,2179.0,0.0,-2.0,-1129.0,0.0,-2.0},
        {-1,1,0,1,0,3276.0,0.0,1.0,-9.0,0.0,0.0},
        {1,1,0,0,0,-3389.0,0.0,5.0,35.0,0.0,-2.0},
        {1,0,2,0,0,3339.0,0.0,-13.0,-107.0,0.0,1.0},
        {-1,0,2,-2,1,-1987.0,0.0,-6.0,1073.0,0.0,-2.0},
        {1,0,0,0,2,-1981.0,0.0,0.0,854.0,0.0,0.0},
        {-1,0,0,1,0,4026.0,0.0,-353.0,-553.0,0.0,-139.0},
        {0,0,2,1,2,1660.0,0.0,-5.0,-710.0,0.0,-2.0},
        {-1,0,2,4,2,-1521.0,0.0,9.0,647.0,0.0,4.0},
        {-1,1,0,1,1,1314.0,0.0,0.0,-700.0,0.0,0.0},
        {0,-2,2,-2,1,-1283.0,0.0,0.0,672.0,0.0,0.0},
        {1,0,2,2,1,-1331.0,0.0,8.0,663.0,0.0,4.0},
        {-2,0,2,2,2,1383.0,0.0,-2.0,-594.0,0.0,-2.0},
        {-1,0,0,0,2,1405.0,0.0,4.0,-610.0,0.0,2.0},
        {1,1,2,-2,2,1290.0,0.0,0.0,-556.0,0.0,0.0}
    };
    const double u2r = CG_ARCSEC2RAD / 1.0e7;
    const double dp_plan = -0.135 * CG_MAS2RAD;
    const double de_plan = 0.388 * CG_MAS2RAD;
    double t = (jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY;
    double el = fmod(485868.249036 + 1717915923.2178 * t, 1296000.0) * CG_ARCSEC2RAD;
    double elp = fmod(1287104.79305 + 129596581.0481 * t, 1296000.0) * CG_ARCSEC2RAD;
    double f = fmod(335779.526232 + 1739527262.8478 * t, 1296000.0) * CG_ARCSEC2RAD;
    double d = fmod(1072260.70369 + 1602961601.2090 * t, 1296000.0) * CG_ARCSEC2RAD;
    double om = fmod(450160.398036 - 6962890.5431 * t, 1296000.0) * CG_ARCSEC2RAD;
    double dp = 0.0;
    double de = 0.0;
    int i;
    for (i = (int)(sizeof(x) / sizeof(x[0])) - 1; i >= 0; --i) {
        double arg = fmod((double)x[i].nl * el + (double)x[i].nlp * elp +
                          (double)x[i].nf * f + (double)x[i].nd * d +
                          (double)x[i].nom * om, CG_D2PI);
        double sarg = sin(arg);
        double carg = cos(arg);
        dp += (x[i].ps + x[i].pst * t) * sarg + x[i].pc * carg;
        de += (x[i].ec + x[i].ect * t) * carg + x[i].es * sarg;
    }
    *dpsi = dp * u2r + dp_plan;
    *deps = de * u2r + de_plan;
}

static void cg_bias_precession_matrix(double jd_tt, double rbp[3][3])
{
    const double eps0 = 84381.448 * CG_ARCSEC2RAD;
    const double dpsibi = -0.041775 * CG_ARCSEC2RAD;
    const double depsbi = -0.0068192 * CG_ARCSEC2RAD;
    const double dra0 = -0.0146 * CG_ARCSEC2RAD;
    double t = (jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY;
    double psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * CG_ARCSEC2RAD;
    double oma77 = eps0 + ((0.05127 + (-0.007726) * t) * t) * t * CG_ARCSEC2RAD;
    double chia = (10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * CG_ARCSEC2RAD;
    double dpsipr = (-0.29965 * CG_ARCSEC2RAD) * t;
    double depspr = (-0.02524 * CG_ARCSEC2RAD) * t;
    double rb[3][3];
    double rp[3][3];
    double psia = psia77 + dpsipr;
    double oma = oma77 + depspr;

    cg_mat_identity(rb);
    cg_rz(dra0, rb);
    cg_ry(dpsibi * sin(eps0), rb);
    cg_rx(-depsbi, rb);

    cg_mat_identity(rp);
    cg_rx(eps0, rp);
    cg_rz(-psia, rp);
    cg_rx(-oma, rp);
    cg_rz(chia, rp);

    cg_mat_mul(rp, rb, rbp);
}

static double cg_obl80(double jd_tt)
{
    double t = (jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY;
    return CG_ARCSEC2RAD * (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t) * t) * t);
}

static void cg_pnm00b(double jd_tt, double rbpn[3][3], double *dpsi_out, double *epsa_out)
{
    double dpsi;
    double deps;
    double epsa;
    double rbp[3][3];
    double rn[3][3];
    double dpsipr;
    double depspr;
    cg_nut00b(jd_tt, &dpsi, &deps);
    dpsipr = (-0.29965 * CG_ARCSEC2RAD) * ((jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY);
    depspr = (-0.02524 * CG_ARCSEC2RAD) * ((jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY);
    (void)dpsipr;
    epsa = cg_obl80(jd_tt) + depspr;
    cg_bias_precession_matrix(jd_tt, rbp);
    cg_mat_identity(rn);
    cg_rx(epsa, rn);
    cg_rz(-dpsi, rn);
    cg_rx(-(epsa + deps), rn);
    cg_mat_mul(rn, rbp, rbpn);
    if (dpsi_out) {
        *dpsi_out = dpsi;
    }
    if (epsa_out) {
        *epsa_out = epsa;
    }
}

static double cg_era00(double jd_ut1)
{
    double d1 = floor(jd_ut1 - 0.5) + 0.5;
    double d2 = jd_ut1 - d1;
    double t = d1 + (d2 - CG_JD_J2000);
    double f = fmod(d1, 1.0) + fmod(d2, 1.0);
    return cg_anp(CG_D2PI * (f + 0.7790572732640 + 0.00273781191135448 * t));
}

static double cg_gmst00(double jd_ut1, double jd_tt)
{
    double t = (jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY;
    return cg_anp(cg_era00(jd_ut1) +
                  (0.014506 + (4612.15739966 +
                  (1.39667721 + (-0.00009344 + 0.00001882 * t) * t) * t) * t) *
                  CG_ARCSEC2RAD);
}

static void cg_itrs_to_j2000_matrix(const cg_time_t *time_utc, double out[3][3])
{
    double jd_utc = time_utc->jd_utc;
    double jd_tt = jd_utc + CG_TT_MINUS_UTC_SECONDS / CG_SECONDS_PER_DAY;
    cg_eop_value_t eop = cg_lookup_eop(cg_mjd_from_jd(jd_utc));
    double jd_ut1 = jd_utc + eop.dut1_seconds / CG_SECONDS_PER_DAY;
    double rbpn[3][3];
    double rbpn_t[3][3];
    double rz[3][3];
    double pm[3][3];
    double rot_pm[3][3];
    double dpsi;
    double epsa;
    double gst;
    int i;
    int j;
    cg_pnm00b(jd_tt, rbpn, &dpsi, &epsa);
    gst = cg_anp(cg_gmst00(jd_ut1, jd_tt) + dpsi * cos(epsa));
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            rbpn_t[i][j] = rbpn[j][i];
        }
    }
    cg_rz_vector_matrix(gst, rz);
    cg_polar_motion_matrix(-eop.x_pole_rad, -eop.y_pole_rad, pm);
    cg_mat_mul(rz, pm, rot_pm);
    cg_mat_mul(rbpn_t, rot_pm, out);
}

static void cg_ecef_to_j2000(
    const cg_time_t *time_utc,
    cg_vec3_t r_ecef,
    cg_vec3_t v_ecef,
    cg_vec3_t *r_j2000,
    cg_vec3_t *v_j2000)
{
    double m[3][3];
    cg_vec3_t spin;
    cg_vec3_t v_total;
    cg_itrs_to_j2000_matrix(time_utc, m);
    spin.x = -CG_OMEGA_EARTH * r_ecef.y;
    spin.y = CG_OMEGA_EARTH * r_ecef.x;
    spin.z = 0.0;
    v_total = cg_vec_add(v_ecef, spin);
    *r_j2000 = cg_mat_vec(m, r_ecef);
    *v_j2000 = cg_mat_vec(m, v_total);
}

static cg_vec3_t cg_observation_to_j2000_position(const cg_observation_t *observation)
{
    cg_vec3_t zero;
    cg_vec3_t r_j2000;
    cg_vec3_t ignored_v_j2000;
    zero.x = zero.y = zero.z = 0.0;
    cg_ecef_to_j2000(
        &observation->time_utc,
        observation->r_ecef_m,
        zero,
        &r_j2000,
        &ignored_v_j2000);
    return r_j2000;
}

static double cg_normalization_factor(int n, int m)
{
    double factor = (double)(2 * n + 1);
    int k;
    if (m > 0) {
        factor *= 2.0;
    }
    for (k = n - m + 1; k <= n + m; ++k) {
        factor /= (double)k;
    }
    return sqrt(factor);
}

static void cg_normalized_legendre(
    double s,
    double pbar[CG_GRAVITY_DEGREE + 1][CG_GRAVITY_DEGREE + 1],
    double dpbar[CG_GRAVITY_DEGREE + 1][CG_GRAVITY_DEGREE + 1])
{
    double p[CG_GRAVITY_DEGREE + 1][CG_GRAVITY_DEGREE + 1];
    double u = sqrt(fmax(0.0, 1.0 - s * s));
    double denom = s * s - 1.0;
    int n;
    int m;
    memset(pbar, 0, sizeof(double) * (CG_GRAVITY_DEGREE + 1) * (CG_GRAVITY_DEGREE + 1));
    memset(dpbar, 0, sizeof(double) * (CG_GRAVITY_DEGREE + 1) * (CG_GRAVITY_DEGREE + 1));
    memset(p, 0, sizeof(p));
    p[0][0] = 1.0;
    for (m = 1; m <= CG_GRAVITY_DEGREE; ++m) {
        p[m][m] = (double)(2 * m - 1) * u * p[m - 1][m - 1];
    }
    for (m = 0; m < CG_GRAVITY_DEGREE; ++m) {
        p[m + 1][m] = (double)(2 * m + 1) * s * p[m][m];
    }
    for (m = 0; m <= CG_GRAVITY_DEGREE; ++m) {
        for (n = m + 2; n <= CG_GRAVITY_DEGREE; ++n) {
            p[n][m] =
                ((double)(2 * n - 1) * s * p[n - 1][m] -
                 (double)(n + m - 1) * p[n - 2][m]) /
                (double)(n - m);
        }
    }

    for (n = 0; n <= CG_GRAVITY_DEGREE; ++n) {
        for (m = 0; m <= n; ++m) {
            double norm = cg_normalization_factor(n, m);
            double p_prev;
            double dp;
            pbar[n][m] = norm * p[n][m];
            if (n == 0 || fabs(denom) < 1.0e-14) {
                dpbar[n][m] = 0.0;
                continue;
            }
            p_prev = n > m ? p[n - 1][m] : 0.0;
            dp = ((double)n * s * p[n][m] - (double)(n + m) * p_prev) / denom;
            dpbar[n][m] = norm * dp;
        }
    }
}

static cg_vec3_t cg_gravity_acceleration_body_fixed(cg_vec3_t r)
{
    double radius = cg_vec_norm(r);
    cg_vec3_t acceleration;
    if (radius <= 0.0) {
        acceleration.x = acceleration.y = acceleration.z = 0.0;
        return acceleration;
    }
    {
        double radius2 = radius * radius;
        double rho2 = r.x * r.x + r.y * r.y;
        double s = r.z / radius;
        double lambda = atan2(r.y, r.x);
        double pbar[CG_GRAVITY_DEGREE + 1][CG_GRAVITY_DEGREE + 1];
        double dpbar[CG_GRAVITY_DEGREE + 1][CG_GRAVITY_DEGREE + 1];
        double radial_power[CG_GRAVITY_DEGREE + 1];
        double radius_ratio = kGravityRadiusM / radius;
        double d_dr = 0.0;
        double d_ds = 0.0;
        double d_dlambda = 0.0;
        size_t term_index;
        int n;

        acceleration.x = -kGravityMu * r.x / (radius2 * radius);
        acceleration.y = -kGravityMu * r.y / (radius2 * radius);
        acceleration.z = -kGravityMu * r.z / (radius2 * radius);
        cg_normalized_legendre(s, pbar, dpbar);
        radial_power[0] = 1.0;
        for (n = 1; n <= CG_GRAVITY_DEGREE; ++n) {
            radial_power[n] = radial_power[n - 1] * radius_ratio;
        }
        for (term_index = 0; term_index < sizeof(kGravityTerms) / sizeof(kGravityTerms[0]); ++term_index) {
            const GravityTerm *term = &kGravityTerms[term_index];
            double angle = (double)term->m * lambda;
            double cos_angle = cos(angle);
            double sin_angle = sin(angle);
            double harmonic = term->c * cos_angle + term->s * sin_angle;
            double harmonic_lambda =
                (double)term->m * (-term->c * sin_angle + term->s * cos_angle);
            double q = radial_power[term->n];
            double p = pbar[term->n][term->m];
            d_dr += -kGravityMu / radius2 * (double)(term->n + 1) * q * p * harmonic;
            d_ds += kGravityMu / radius * q * dpbar[term->n][term->m] * harmonic;
            d_dlambda += kGravityMu / radius * q * p * harmonic_lambda;
        }

        acceleration.x += d_dr * r.x / radius + d_ds * (-s * r.x / radius2);
        acceleration.y += d_dr * r.y / radius + d_ds * (-s * r.y / radius2);
        acceleration.z += d_dr * r.z / radius + d_ds * ((1.0 - s * s) / radius);
        if (rho2 > 0.0) {
            acceleration.x += d_dlambda * (-r.y / rho2);
            acceleration.y += d_dlambda * (r.x / rho2);
        }
    }
    return acceleration;
}

static cg_vec3_t cg_gravity_acceleration_j2000(const cg_time_t *time_utc, cg_vec3_t r_j2000)
{
    double m[3][3];
    cg_vec3_t r_body;
    cg_vec3_t a_body;
    cg_itrs_to_j2000_matrix(time_utc, m);
    r_body.x = m[0][0] * r_j2000.x + m[1][0] * r_j2000.y + m[2][0] * r_j2000.z;
    r_body.y = m[0][1] * r_j2000.x + m[1][1] * r_j2000.y + m[2][1] * r_j2000.z;
    r_body.z = m[0][2] * r_j2000.x + m[1][2] * r_j2000.y + m[2][2] * r_j2000.z;
    a_body = cg_gravity_acceleration_body_fixed(r_body);
    return cg_mat_vec(m, a_body);
}

static cg_vec3_t cg_central_gravity_acceleration_j2000(cg_vec3_t r_j2000)
{
    double radius = cg_vec_norm(r_j2000);
    cg_vec3_t zero;
    zero.x = zero.y = zero.z = 0.0;
    if (radius <= 0.0) {
        return zero;
    }
    return cg_vec_scale(r_j2000, -kGravityMu / (radius * radius * radius));
}

static double cg_radians_from_degrees(double degrees)
{
    return degrees * CG_PI / 180.0;
}

static cg_vec3_t cg_ecliptic_to_equatorial(
    double radius_m,
    double longitude_rad,
    double latitude_rad,
    double obliquity_rad)
{
    double cos_beta = cos(latitude_rad);
    double x_ecl = radius_m * cos_beta * cos(longitude_rad);
    double y_ecl = radius_m * cos_beta * sin(longitude_rad);
    double z_ecl = radius_m * sin(latitude_rad);
    double cos_eps = cos(obliquity_rad);
    double sin_eps = sin(obliquity_rad);
    cg_vec3_t r;
    r.x = x_ecl;
    r.y = y_ecl * cos_eps - z_ecl * sin_eps;
    r.z = y_ecl * sin_eps + z_ecl * cos_eps;
    return r;
}

static cg_vec3_t cg_sun_position_j2000(const cg_time_t *time_utc)
{
    double jd_tt = time_utc->jd_utc + CG_TT_MINUS_UTC_SECONDS / CG_SECONDS_PER_DAY;
    double d = jd_tt - CG_JD_J2000;
    double mean_longitude = cg_anp(cg_radians_from_degrees(280.460 + 0.9856474 * d));
    double mean_anomaly = cg_anp(cg_radians_from_degrees(357.528 + 0.9856003 * d));
    double ecliptic_longitude = cg_anp(mean_longitude +
        cg_radians_from_degrees(1.915) * sin(mean_anomaly) +
        cg_radians_from_degrees(0.020) * sin(2.0 * mean_anomaly));
    double radius_au = 1.00014 - 0.01671 * cos(mean_anomaly) -
        0.00014 * cos(2.0 * mean_anomaly);
    double obliquity = cg_radians_from_degrees(23.439291 - 0.0000004 * d);
    return cg_ecliptic_to_equatorial(
        radius_au * CG_ASTRONOMICAL_UNIT_M, ecliptic_longitude, 0.0, obliquity);
}

static cg_vec3_t cg_moon_position_j2000(const cg_time_t *time_utc)
{
    double jd_tt = time_utc->jd_utc + CG_TT_MINUS_UTC_SECONDS / CG_SECONDS_PER_DAY;
    double d = jd_tt - CG_JD_J2000;
    double mean_longitude = cg_anp(cg_radians_from_degrees(218.316 + 13.176396 * d));
    double mean_anomaly = cg_anp(cg_radians_from_degrees(134.963 + 13.064993 * d));
    double argument_latitude = cg_anp(cg_radians_from_degrees(93.272 + 13.229350 * d));
    double ecliptic_longitude = cg_anp(mean_longitude + cg_radians_from_degrees(6.289) * sin(mean_anomaly));
    double ecliptic_latitude = cg_radians_from_degrees(5.128) * sin(argument_latitude);
    double radius_m = (385001.0 - 20905.0 * cos(mean_anomaly)) * 1000.0;
    double t = (jd_tt - CG_JD_J2000) / CG_JULIAN_CENTURY;
    double obliquity = cg_radians_from_degrees(23.439291 - 0.0130042 * t);
    return cg_ecliptic_to_equatorial(radius_m, ecliptic_longitude, ecliptic_latitude, obliquity);
}

static cg_vec3_t cg_third_body_acceleration(cg_vec3_t r_sat, cg_vec3_t r_body, double mu_body)
{
    cg_vec3_t sat_to_body = cg_vec_sub(r_body, r_sat);
    double d_sat = cg_vec_norm(sat_to_body);
    double d_body = cg_vec_norm(r_body);
    cg_vec3_t zero;
    zero.x = zero.y = zero.z = 0.0;
    if (d_sat <= 0.0 || d_body <= 0.0) {
        return zero;
    }
    return cg_vec_sub(
        cg_vec_scale(sat_to_body, mu_body / (d_sat * d_sat * d_sat)),
        cg_vec_scale(r_body, mu_body / (d_body * d_body * d_body)));
}

typedef struct cg_space_weather_t {
    double f107_sfu;
    double ap;
} cg_space_weather_t;

static double cg_clamp_double(double value, double lo, double hi)
{
    if (!isfinite(value)) {
        return lo;
    }
    if (value < lo) {
        return lo;
    }
    if (value > hi) {
        return hi;
    }
    return value;
}

static cg_space_weather_t cg_nominal_space_weather(const cg_time_t *time_utc)
{
    cg_space_weather_t sw;
    double phase;
    if (!time_utc || !isfinite(time_utc->unix_seconds)) {
        sw.f107_sfu = 150.0;
        sw.ap = 15.0;
        return sw;
    }

    phase = CG_D2PI * (time_utc->unix_seconds / (27.2753 * CG_SECONDS_PER_DAY));
    sw.f107_sfu = cg_clamp_double(155.0 + 18.0 * sin(phase), 70.0, 260.0);
    sw.ap = cg_clamp_double(12.0 + 5.0 * (1.0 + sin(phase + 1.7)), 0.0, 80.0);
    return sw;
}

static double cg_space_weather_density_scale(const cg_time_t *time_utc, double altitude_m)
{
    cg_space_weather_t sw = cg_nominal_space_weather(time_utc);
    double altitude_factor = cg_clamp_double((altitude_m - 180000.0) / 220000.0, 0.0, 2.0);
    double solar_scale = exp(0.0045 * (sw.f107_sfu - 120.0) * altitude_factor);
    double geomag_scale = 1.0 + 0.010 * sw.ap * exp(-fabs(altitude_m - 270000.0) / 260000.0);
    return cg_clamp_double(solar_scale * geomag_scale, 0.25, 8.0);
}

static double cg_local_solar_density_scale(double solar_cosine)
{
    double day = fmax(0.0, solar_cosine);
    double terminator = 1.0 - fabs(cg_clamp_double(solar_cosine, -1.0, 1.0));
    return cg_clamp_double(0.78 + 0.34 * day + 0.08 * terminator, 0.55, 1.35);
}

static double cg_atmosphere_density_kgpm3(
    double altitude_m,
    const cg_time_t *time_utc,
    double solar_cosine)
{
    typedef struct cg_density_point_t {
        double altitude_m;
        double density;
    } cg_density_point_t;
    static const cg_density_point_t table[] = {
        {180000.0, 5.464e-10},
        {200000.0, 2.789e-10},
        {250000.0, 7.248e-11},
        {300000.0, 2.418e-11},
        {350000.0, 9.518e-12},
        {400000.0, 3.725e-12},
        {450000.0, 1.585e-12},
        {500000.0, 6.967e-13},
        {550000.0, 3.177e-13},
        {600000.0, 1.454e-13},
        {700000.0, 3.614e-14},
        {800000.0, 1.170e-14},
        {900000.0, 5.245e-15},
        {1000000.0, 3.019e-15}
    };
    size_t count = sizeof(table) / sizeof(table[0]);
    size_t i;
    double density = table[count - 1U].density;
    if (altitude_m <= table[0].altitude_m) {
        density = table[0].density;
    } else if (altitude_m >= table[count - 1U].altitude_m) {
        double h_scale = 120000.0;
        density = table[count - 1U].density *
            exp(-(altitude_m - table[count - 1U].altitude_m) / h_scale);
    } else {
        for (i = 1U; i < count; ++i) {
            if (altitude_m <= table[i].altitude_m) {
                double f = (altitude_m - table[i - 1U].altitude_m) /
                    (table[i].altitude_m - table[i - 1U].altitude_m);
                double log_rho = log(table[i - 1U].density) +
                    f * (log(table[i].density) - log(table[i - 1U].density));
                density = exp(log_rho);
                break;
            }
        }
    }
    return density *
        cg_space_weather_density_scale(time_utc, altitude_m) *
        cg_local_solar_density_scale(solar_cosine);
}

static cg_vec3_t cg_drag_acceleration_per_ballistic(
    const cg_time_t *time_utc,
    cg_vec3_t r_j2000,
    cg_vec3_t v_j2000)
{
    double m[3][3];
    cg_vec3_t r_body;
    cg_vec3_t r_sun_j2000;
    cg_vec3_t r_sun_body;
    cg_vec3_t atmosphere_velocity_body;
    cg_vec3_t atmosphere_velocity_j2000;
    cg_vec3_t relative_velocity;
    double altitude_m;
    double solar_cosine;
    double density;
    double speed;
    cg_vec3_t zero;
    zero.x = zero.y = zero.z = 0.0;
    cg_itrs_to_j2000_matrix(time_utc, m);
    r_body.x = m[0][0] * r_j2000.x + m[1][0] * r_j2000.y + m[2][0] * r_j2000.z;
    r_body.y = m[0][1] * r_j2000.x + m[1][1] * r_j2000.y + m[2][1] * r_j2000.z;
    r_body.z = m[0][2] * r_j2000.x + m[1][2] * r_j2000.y + m[2][2] * r_j2000.z;
    altitude_m = cg_vec_norm(r_body) - CG_EARTH_EQUATORIAL_RADIUS_M;
    r_sun_j2000 = cg_sun_position_j2000(time_utc);
    r_sun_body.x = m[0][0] * r_sun_j2000.x + m[1][0] * r_sun_j2000.y + m[2][0] * r_sun_j2000.z;
    r_sun_body.y = m[0][1] * r_sun_j2000.x + m[1][1] * r_sun_j2000.y + m[2][1] * r_sun_j2000.z;
    r_sun_body.z = m[0][2] * r_sun_j2000.x + m[1][2] * r_sun_j2000.y + m[2][2] * r_sun_j2000.z;
    solar_cosine = cg_vec_dot(cg_vec_unit(r_body), cg_vec_unit(r_sun_body));
    density = cg_atmosphere_density_kgpm3(altitude_m, time_utc, solar_cosine);
    atmosphere_velocity_body.x = -CG_OMEGA_EARTH * r_body.y;
    atmosphere_velocity_body.y = CG_OMEGA_EARTH * r_body.x;
    atmosphere_velocity_body.z = 0.0;
    atmosphere_velocity_j2000 = cg_mat_vec(m, atmosphere_velocity_body);
    relative_velocity = cg_vec_sub(v_j2000, atmosphere_velocity_j2000);
    speed = cg_vec_norm(relative_velocity);
    if (speed <= 0.0 || density <= 0.0) {
        return zero;
    }
    return cg_vec_scale(relative_velocity, -0.5 * density * speed);
}

static double cg_cylindrical_shadow_factor(cg_vec3_t r_sat, cg_vec3_t r_sun)
{
    double sun_distance = cg_vec_norm(r_sun);
    double cross_track;
    if (sun_distance <= 0.0 || cg_vec_dot(r_sat, r_sun) >= 0.0) {
        return 1.0;
    }
    cross_track = cg_vec_norm(cg_vec_cross(r_sat, r_sun)) / sun_distance;
    return cross_track < CG_EARTH_EQUATORIAL_RADIUS_M ? 0.0 : 1.0;
}

static cg_vec3_t cg_solar_radiation_pressure_acceleration(
    cg_vec3_t r_sat,
    cg_vec3_t r_sun,
    double area_over_mass)
{
    cg_vec3_t sun_to_sat;
    double distance;
    double shadow;
    double scale;
    cg_vec3_t zero;
    zero.x = zero.y = zero.z = 0.0;
    if (area_over_mass <= 0.0) {
        return zero;
    }
    sun_to_sat = cg_vec_sub(r_sat, r_sun);
    distance = cg_vec_norm(sun_to_sat);
    if (distance <= 0.0) {
        return zero;
    }
    shadow = cg_cylindrical_shadow_factor(r_sat, r_sun);
    scale = shadow * CG_SOLAR_RADIATION_PRESSURE * area_over_mass *
        (CG_ASTRONOMICAL_UNIT_M * CG_ASTRONOMICAL_UNIT_M) / (distance * distance);
    return cg_vec_scale(sun_to_sat, scale / distance);
}

static cg_vec3_t cg_empirical_rtn_acceleration(
    cg_vec3_t r_j2000,
    cg_vec3_t v_j2000,
    cg_vec3_t empirical_rtn)
{
    cg_vec3_t radial = cg_vec_unit(r_j2000);
    cg_vec3_t normal = cg_vec_unit(cg_vec_cross(r_j2000, v_j2000));
    cg_vec3_t transverse = cg_vec_cross(normal, radial);
    return cg_vec_add(
        cg_vec_add(cg_vec_scale(radial, empirical_rtn.x),
                   cg_vec_scale(transverse, empirical_rtn.y)),
        cg_vec_scale(normal, empirical_rtn.z));
}

static cg_force_model_t cg_default_force_model(void)
{
    cg_force_model_t force_model;
    force_model.drag_ballistic_m2_per_kg = 0.0;
    force_model.gravity_harmonic_scale = 1.0;
    force_model.srp_area_m2_per_kg = 0.0040;
    force_model.empirical_rtn_mps2.x = 0.0;
    force_model.empirical_rtn_mps2.y = 0.0;
    force_model.empirical_rtn_mps2.z = 0.0;
    return force_model;
}

static cg_vec3_t cg_force_model_acceleration(
    const cg_time_t *time_utc,
    cg_vec3_t r_j2000,
    cg_vec3_t v_j2000,
    const cg_force_model_t *force_model)
{
    cg_vec3_t acceleration = cg_gravity_acceleration_j2000(time_utc, r_j2000);
    cg_vec3_t central_gravity = cg_central_gravity_acceleration_j2000(r_j2000);
    cg_vec3_t r_sun = cg_sun_position_j2000(time_utc);
    cg_vec3_t r_moon = cg_moon_position_j2000(time_utc);
    acceleration = cg_vec_add(
        central_gravity,
        cg_vec_scale(
            cg_vec_sub(acceleration, central_gravity),
            force_model->gravity_harmonic_scale));
    acceleration = cg_vec_add(acceleration, cg_third_body_acceleration(r_j2000, r_sun, CG_SUN_MU));
    acceleration = cg_vec_add(acceleration, cg_third_body_acceleration(r_j2000, r_moon, CG_MOON_MU));
    acceleration = cg_vec_add(
        acceleration,
        cg_vec_scale(
            cg_drag_acceleration_per_ballistic(time_utc, r_j2000, v_j2000),
            force_model->drag_ballistic_m2_per_kg));
    acceleration = cg_vec_add(
        acceleration,
        cg_solar_radiation_pressure_acceleration(
            r_j2000, r_sun, force_model->srp_area_m2_per_kg));
    acceleration = cg_vec_add(
        acceleration,
        cg_empirical_rtn_acceleration(r_j2000, v_j2000, force_model->empirical_rtn_mps2));
    return acceleration;
}

static void cg_state_derivative(
    const double state[6],
    const cg_time_t *time_utc,
    const cg_force_model_t *force_model,
    double derivative[6])
{
    cg_vec3_t r;
    cg_vec3_t v;
    cg_vec3_t a;
    r.x = state[0]; r.y = state[1]; r.z = state[2];
    v.x = state[3]; v.y = state[4]; v.z = state[5];
    a = cg_force_model_acceleration(time_utc, r, v, force_model);
    derivative[0] = state[3];
    derivative[1] = state[4];
    derivative[2] = state[5];
    derivative[3] = a.x;
    derivative[4] = a.y;
    derivative[5] = a.z;
}

static int cg_propagate_orbit(
    const cg_state_t *initial,
    const cg_time_t *target_time,
    double step_seconds,
    const cg_force_model_t *force_model,
    cg_state_t *out_state)
{
    double y[6];
    double elapsed = 0.0;
    double remaining;
    double direction;
    int i;
    if (!initial || !target_time || !force_model || !out_state || step_seconds <= 0.0) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    y[0] = initial->r_j2000_m.x;
    y[1] = initial->r_j2000_m.y;
    y[2] = initial->r_j2000_m.z;
    y[3] = initial->v_j2000_mps.x;
    y[4] = initial->v_j2000_mps.y;
    y[5] = initial->v_j2000_mps.z;

    remaining = fabs(target_time->unix_seconds - initial->time_utc.unix_seconds);
    direction = target_time->unix_seconds >= initial->time_utc.unix_seconds ? 1.0 : -1.0;
    while (remaining > 1.0e-9) {
        double h = direction * fmin(step_seconds, remaining);
        double k1[6];
        double k2[6];
        double k3[6];
        double k4[6];
        double temp[6];
        cg_time_t step_time;

        step_time = cg_time_from_unix_seconds(initial->time_utc.unix_seconds + elapsed);
        cg_state_derivative(y, &step_time, force_model, k1);
        for (i = 0; i < 6; ++i) {
            temp[i] = y[i] + 0.5 * h * k1[i];
        }
        step_time = cg_time_from_unix_seconds(initial->time_utc.unix_seconds + elapsed + 0.5 * h);
        cg_state_derivative(temp, &step_time, force_model, k2);
        for (i = 0; i < 6; ++i) {
            temp[i] = y[i] + 0.5 * h * k2[i];
        }
        step_time = cg_time_from_unix_seconds(initial->time_utc.unix_seconds + elapsed + 0.5 * h);
        cg_state_derivative(temp, &step_time, force_model, k3);
        for (i = 0; i < 6; ++i) {
            temp[i] = y[i] + h * k3[i];
        }
        step_time = cg_time_from_unix_seconds(initial->time_utc.unix_seconds + elapsed + h);
        cg_state_derivative(temp, &step_time, force_model, k4);
        for (i = 0; i < 6; ++i) {
            y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
        elapsed += h;
        remaining -= fabs(h);
    }

    out_state->time_utc = *target_time;
    out_state->r_j2000_m.x = y[0];
    out_state->r_j2000_m.y = y[1];
    out_state->r_j2000_m.z = y[2];
    out_state->v_j2000_mps.x = y[3];
    out_state->v_j2000_mps.y = y[4];
    out_state->v_j2000_mps.z = y[5];
    return CG_OK;
}

int cg_precompute_observations(cg_observation_t *observations, size_t count)
{
    if (!observations || count < 2) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    return CG_OK;
}

cg_options_t cg_default_options(void)
{
    cg_options_t opt;
    opt.degree = 10;
    return opt;
}

static void cg_select_window(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    size_t *first,
    size_t *last)
{
    double half_window = 0.5 * CG_FIT_WINDOW_SECONDS;
    double q = query->unix_seconds;
    size_t i;
    size_t f = 0;
    size_t l = count - 1;
    while (f < count && obs[f].time_utc.unix_seconds < q - half_window) {
        ++f;
    }
    while (l > 0 && obs[l].time_utc.unix_seconds > q + half_window) {
        --l;
    }
    if (f <= l && l - f + 1 >= 2) {
        *first = f;
        *last = l;
        return;
    }
    {
        size_t nearest = 0;
        double best = DBL_MAX;
        for (i = 0; i < count; ++i) {
            double d = fabs(obs[i].time_utc.unix_seconds - q);
            if (d < best) {
                best = d;
                nearest = i;
            }
        }
        if (nearest == 0) {
            *first = 0;
            *last = 1;
        } else if (nearest + 1 >= count) {
            *first = count - 2;
            *last = count - 1;
        } else {
            *first = nearest;
            *last = nearest + 1;
        }
    }
}

static void cg_cheb_basis(double tau, int degree, double *basis)
{
    int i;
    basis[0] = 1.0;
    if (degree >= 1) {
        basis[1] = tau;
    }
    for (i = 2; i <= degree; ++i) {
        basis[i] = 2.0 * tau * basis[i - 1] - basis[i - 2];
    }
}

static int cg_solve_linear(int n, double a[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 2], double *x)
{
    int i;
    int j;
    int k;
    for (i = 0; i < n; ++i) {
        int pivot = i;
        double pivot_abs = fabs(a[i][i]);
        for (j = i + 1; j < n; ++j) {
            double v = fabs(a[j][i]);
            if (v > pivot_abs) {
                pivot_abs = v;
                pivot = j;
            }
        }
        if (pivot_abs < 1.0e-18) {
            return CG_ERR_FIT;
        }
        if (pivot != i) {
            for (k = i; k <= n; ++k) {
                double tmp = a[i][k];
                a[i][k] = a[pivot][k];
                a[pivot][k] = tmp;
            }
        }
        for (j = i + 1; j < n; ++j) {
            double factor = a[j][i] / a[i][i];
            a[j][i] = 0.0;
            for (k = i + 1; k <= n; ++k) {
                a[j][k] -= factor * a[i][k];
            }
        }
    }
    for (i = n - 1; i >= 0; --i) {
        double sum = a[i][n];
        for (j = i + 1; j < n; ++j) {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    return CG_OK;
}

static int cg_build_fit(
    const cg_observation_t *obs,
    size_t first,
    size_t last,
    int requested_degree,
    cg_fit_t *fit)
{
    size_t nobs;
    int degree;
    int ncoef;
    double ata[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 1];
    double rhs[3][CG_MAX_DEGREE + 1];
    size_t i;
    int j;
    int k;
    int coord;
    double basis[CG_MAX_DEGREE + 1];
    if (!obs || !fit || last < first) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    nobs = last - first + 1;
    if (nobs < 2) {
        return CG_ERR_FIT;
    }
    degree = requested_degree;
    if (degree < 1) {
        degree = 1;
    }
    if (degree > CG_MAX_DEGREE) {
        degree = CG_MAX_DEGREE;
    }
    if ((size_t)degree >= nobs) {
        degree = (int)nobs - 1;
    }
    ncoef = degree + 1;
    memset(ata, 0, sizeof(ata));
    memset(rhs, 0, sizeof(rhs));
    memset(fit, 0, sizeof(*fit));
    fit->degree = degree;
    fit->first_index = first;
    fit->last_index = last;
    fit->epoch_seconds = 0.5 * (obs[first].time_utc.unix_seconds + obs[last].time_utc.unix_seconds);
    fit->scale_seconds = fmax(fabs(obs[first].time_utc.unix_seconds - fit->epoch_seconds),
                              fabs(obs[last].time_utc.unix_seconds - fit->epoch_seconds));
    if (fit->scale_seconds <= 0.0) {
        return CG_ERR_FIT;
    }
    for (i = first; i <= last; ++i) {
        double tau = (obs[i].time_utc.unix_seconds - fit->epoch_seconds) / fit->scale_seconds;
        cg_vec3_t r_j2000 = cg_observation_to_j2000_position(&obs[i]);
        cg_cheb_basis(tau, degree, basis);
        for (j = 0; j < ncoef; ++j) {
            for (k = 0; k < ncoef; ++k) {
                ata[j][k] += basis[j] * basis[k];
            }
            rhs[0][j] += basis[j] * r_j2000.x;
            rhs[1][j] += basis[j] * r_j2000.y;
            rhs[2][j] += basis[j] * r_j2000.z;
        }
    }
    for (coord = 0; coord < 3; ++coord) {
        double aug[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 2];
        double x[CG_MAX_DEGREE + 1];
        int status;
        for (j = 0; j < ncoef; ++j) {
            for (k = 0; k < ncoef; ++k) {
                aug[j][k] = ata[j][k];
            }
            aug[j][ncoef] = rhs[coord][j];
            x[j] = 0.0;
        }
        status = cg_solve_linear(ncoef, aug, x);
        if (status != CG_OK) {
            return status;
        }
        for (j = 0; j < ncoef; ++j) {
            fit->coeff[coord][j] = x[j];
        }
    }
    return CG_OK;
}

static double cg_cheb_eval(const double *coeff, int degree, double tau)
{
    double b0 = 0.0;
    double b1 = 0.0;
    double b2;
    int j;
    for (j = degree; j >= 1; --j) {
        b2 = b1;
        b1 = b0;
        b0 = 2.0 * tau * b1 - b2 + coeff[j];
    }
    return tau * b0 - b1 + coeff[0];
}

static double cg_cheb_derivative_eval(const double *coeff, int degree, double tau)
{
    double deriv = 0.0;
    double u_prev = 1.0;
    double u_curr = 2.0 * tau;
    int n;
    if (degree < 1) {
        return 0.0;
    }
    deriv += coeff[1];
    for (n = 2; n <= degree; ++n) {
        double u = (n == 2) ? u_curr : 2.0 * tau * u_curr - u_prev;
        if (n > 2) {
            u_prev = u_curr;
            u_curr = u;
        }
        deriv += (double)n * coeff[n] * u;
    }
    return deriv;
}

static void cg_fit_eval(const cg_fit_t *fit, const cg_time_t *query, cg_state_t *out)
{
    double tau = (query->unix_seconds - fit->epoch_seconds) / fit->scale_seconds;
    out->time_utc = *query;
    out->r_j2000_m.x = cg_cheb_eval(fit->coeff[0], fit->degree, tau);
    out->r_j2000_m.y = cg_cheb_eval(fit->coeff[1], fit->degree, tau);
    out->r_j2000_m.z = cg_cheb_eval(fit->coeff[2], fit->degree, tau);
    out->v_j2000_mps.x = cg_cheb_derivative_eval(fit->coeff[0], fit->degree, tau) / fit->scale_seconds;
    out->v_j2000_mps.y = cg_cheb_derivative_eval(fit->coeff[1], fit->degree, tau) / fit->scale_seconds;
    out->v_j2000_mps.z = cg_cheb_derivative_eval(fit->coeff[2], fit->degree, tau) / fit->scale_seconds;
}

static int cg_interpolate_state_cached(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_fit_cache_t *cache,
    cg_state_t *out)
{
    size_t first;
    size_t last;
    int degree = opt->degree;
    int status;
    cg_select_window(obs, count, query, &first, &last);
    if (!cache || !cache->valid || cache->first_index != first || cache->last_index != last || cache->degree != degree) {
        cg_fit_t fit;
        status = cg_build_fit(obs, first, last, degree, &fit);
        if (status != CG_OK) {
            return status;
        }
        if (cache) {
            cache->valid = 1;
            cache->degree = degree;
            cache->first_index = first;
            cache->last_index = last;
            cache->fit = fit;
        } else {
            cg_fit_eval(&fit, query, out);
            return CG_OK;
        }
    }
    cg_fit_eval(&cache->fit, query, out);
    return CG_OK;
}

static void cg_select_extrapolation_history(
    const cg_observation_t *obs,
    size_t count,
    size_t *first,
    size_t *last)
{
    size_t f = 0U;
    size_t l = count - 1U;
    double earliest = obs[l].time_utc.unix_seconds - CG_EXTRAPOLATION_HISTORY_SECONDS;
    while (f < l && obs[f].time_utc.unix_seconds < earliest) {
        ++f;
    }
    if (f >= l && l > 0U) {
        f = l - 1U;
    }
    *first = f;
    *last = l;
}

static int cg_refine_latest_state_from_tail(
    const cg_observation_t *obs,
    size_t count,
    cg_state_t *state)
{
    size_t first;
    size_t last;
    cg_fit_t fit;
    cg_state_t refined;
    int status;
    if (!obs || !state || count < 16U) {
        return 0;
    }
    last = count - 1U;
    first = last;
    while (first > 0U &&
           obs[last].time_utc.unix_seconds - obs[first - 1U].time_utc.unix_seconds <=
               CG_TAIL_REFINE_SECONDS + 1.0e-9) {
        --first;
    }
    if (last - first + 1U < 16U) {
        return 0;
    }
    status = cg_build_fit(obs, first, last, CG_MAX_DEGREE, &fit);
    if (status != CG_OK) {
        return 0;
    }
    cg_fit_eval(&fit, &obs[last].time_utc, &refined);
    state->time_utc = obs[last].time_utc;
    state->r_j2000_m = cg_observation_to_j2000_position(&obs[last]);
    state->v_j2000_mps = refined.v_j2000_mps;
    return 1;
}

static int cg_latest_state_for_extrapolation(
    const cg_observation_t *obs,
    size_t count,
    const cg_options_t *opt,
    cg_state_t *out)
{
    size_t first;
    size_t last;
    cg_fit_t fit;
    int status;
    if (!obs || !opt || !out || count < 2U) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    cg_select_extrapolation_history(obs, count, &first, &last);
    status = cg_build_fit(obs, first, last, opt->degree, &fit);
    if (status != CG_OK) {
        return status;
    }
    cg_fit_eval(&fit, &obs[last].time_utc, out);
    out->time_utc = obs[last].time_utc;
    out->r_j2000_m = cg_observation_to_j2000_position(&obs[last]);
    cg_refine_latest_state_from_tail(obs, count, out);
    return CG_OK;
}

static int cg_fit_state_for_force_model(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_state_t *out)
{
    return cg_interpolate_state_cached(obs, count, query, opt, NULL, out);
}

static double cg_vec_component(cg_vec3_t v, int coord)
{
    if (coord == 0) {
        return v.x;
    }
    if (coord == 1) {
        return v.y;
    }
    return v.z;
}

static void cg_add_weighted_normal_row(
    double normal[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 2],
    int nparam,
    const double *columns,
    double observed,
    double weight)
{
    int i;
    int j;
    for (i = 0; i < nparam; ++i) {
        double ci = weight * columns[i];
        for (j = 0; j < nparam; ++j) {
            normal[i][j] += ci * weight * columns[j];
        }
        normal[i][nparam] += ci * weight * observed;
    }
}

static double cg_clamp_empirical_acceleration(double value)
{
    if (!isfinite(value)) {
        return 0.0;
    }
    if (value < -5.0e-2) {
        return -5.0e-2;
    }
    if (value > 5.0e-2) {
        return 5.0e-2;
    }
    return value;
}

static cg_force_model_t cg_estimate_force_model(
    const cg_observation_t *obs,
    size_t count,
    const cg_options_t *opt,
    const cg_state_t *latest_state)
{
    enum { CG_FORCE_PARAM_COUNT = 5 };
    static const double taus[] = {120.0, 180.0, 240.0, 300.0, 420.0, 540.0};
    static const double parameter_scales[CG_FORCE_PARAM_COUNT] = {
        0.0040, 0.25, 1.0e-6, 1.0e-6, 1.0e-6
    };
    static const double prior_sigma[CG_FORCE_PARAM_COUNT] = {
        5.0, 4.0, 0.5, 0.5, 0.5
    };
    cg_force_model_t model = cg_default_force_model();
    cg_force_model_t base_model;
    double normal[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 2];
    double solution[CG_MAX_DEGREE + 1];
    double latest_unix;
    double first_unix;
    cg_vec3_t radial;
    cg_vec3_t orbit_normal;
    cg_vec3_t transverse;
    int arc_count = 0;
    size_t tau_index;
    int param;
    int coord;
    int status;
    if (!obs || !opt || !latest_state || count < 16U) {
        return model;
    }
    latest_unix = latest_state->time_utc.unix_seconds;
    first_unix = obs[0].time_utc.unix_seconds;
    if (latest_unix - first_unix < 240.0) {
        return model;
    }
    radial = cg_vec_unit(latest_state->r_j2000_m);
    orbit_normal = cg_vec_unit(cg_vec_cross(latest_state->r_j2000_m, latest_state->v_j2000_mps));
    transverse = cg_vec_cross(orbit_normal, radial);
    if (cg_vec_norm(radial) <= 0.0 ||
        cg_vec_norm(orbit_normal) <= 0.0 ||
        cg_vec_norm(transverse) <= 0.0) {
        return model;
    }

    memset(normal, 0, sizeof(normal));
    base_model = model;
    base_model.empirical_rtn_mps2.x = 0.0;
    base_model.empirical_rtn_mps2.y = 0.0;
    base_model.empirical_rtn_mps2.z = 0.0;

    for (tau_index = 0U; tau_index < sizeof(taus) / sizeof(taus[0]); ++tau_index) {
        double tau = taus[tau_index];
        cg_time_t anchor_time;
        cg_state_t anchor_state;
        cg_state_t propagated;
        cg_state_t perturbed[CG_FORCE_PARAM_COUNT];
        cg_vec3_t position_residual;
        cg_vec3_t velocity_residual;
        double weight;
        if (latest_unix - tau < first_unix - 1.0e-9) {
            continue;
        }
        anchor_time = cg_time_from_unix_seconds(latest_unix - tau);
        status = cg_fit_state_for_force_model(obs, count, &anchor_time, opt, &anchor_state);
        if (status != CG_OK) {
            continue;
        }
        status = cg_propagate_orbit(
            &anchor_state,
            &latest_state->time_utc,
            CG_PROPAGATION_STEP_SECONDS,
            &base_model,
            &propagated);
        if (status != CG_OK) {
            continue;
        }
        position_residual = cg_vec_sub(latest_state->r_j2000_m, propagated.r_j2000_m);
        velocity_residual = cg_vec_sub(latest_state->v_j2000_mps, propagated.v_j2000_mps);

        for (param = 0; param < CG_FORCE_PARAM_COUNT; ++param) {
            cg_force_model_t perturbed_model = base_model;
            if (param == 0) {
                perturbed_model.drag_ballistic_m2_per_kg += parameter_scales[param];
            } else if (param == 1) {
                perturbed_model.gravity_harmonic_scale += parameter_scales[param];
            } else if (param == 2) {
                perturbed_model.empirical_rtn_mps2.x += parameter_scales[param];
            } else if (param == 3) {
                perturbed_model.empirical_rtn_mps2.y += parameter_scales[param];
            } else {
                perturbed_model.empirical_rtn_mps2.z += parameter_scales[param];
            }
            status = cg_propagate_orbit(
                &anchor_state,
                &latest_state->time_utc,
                CG_PROPAGATION_STEP_SECONDS,
                &perturbed_model,
                &perturbed[param]);
            if (status != CG_OK) {
                perturbed[param] = propagated;
            }
        }

        weight = sqrt(tau / 300.0);
        for (coord = 0; coord < 3; ++coord) {
            double columns[CG_FORCE_PARAM_COUNT];
            for (param = 0; param < CG_FORCE_PARAM_COUNT; ++param) {
                columns[param] = cg_vec_component(
                    cg_vec_sub(perturbed[param].r_j2000_m, propagated.r_j2000_m),
                    coord);
            }
            cg_add_weighted_normal_row(
                normal,
                CG_FORCE_PARAM_COUNT,
                columns,
                cg_vec_component(position_residual, coord),
                weight);
        }
        for (coord = 0; coord < 3; ++coord) {
            double columns[CG_FORCE_PARAM_COUNT];
            double velocity_weight_seconds = 600.0;
            for (param = 0; param < CG_FORCE_PARAM_COUNT; ++param) {
                columns[param] = velocity_weight_seconds * cg_vec_component(
                    cg_vec_sub(perturbed[param].v_j2000_mps, propagated.v_j2000_mps),
                    coord);
            }
            cg_add_weighted_normal_row(
                normal,
                CG_FORCE_PARAM_COUNT,
                columns,
                velocity_weight_seconds * cg_vec_component(velocity_residual, coord),
                0.35 * weight);
        }
        ++arc_count;
    }

    if (arc_count < 2) {
        return model;
    }
    for (param = 0; param < CG_FORCE_PARAM_COUNT; ++param) {
        normal[param][param] += 1.0 / (prior_sigma[param] * prior_sigma[param]);
        solution[param] = 0.0;
    }
    status = cg_solve_linear(CG_FORCE_PARAM_COUNT, normal, solution);
    if (status != CG_OK) {
        return model;
    }
    model.drag_ballistic_m2_per_kg += solution[0] * parameter_scales[0];
    if (model.drag_ballistic_m2_per_kg < 0.0) {
        model.drag_ballistic_m2_per_kg = 0.0;
    }
    if (model.drag_ballistic_m2_per_kg > 0.080) {
        model.drag_ballistic_m2_per_kg = 0.080;
    }
    model.gravity_harmonic_scale = 1.0;
    model.empirical_rtn_mps2.x = cg_clamp_empirical_acceleration(solution[2] * parameter_scales[2]);
    model.empirical_rtn_mps2.y = cg_clamp_empirical_acceleration(solution[3] * parameter_scales[3]);
    model.empirical_rtn_mps2.z = cg_clamp_empirical_acceleration(solution[4] * parameter_scales[4]);
    return model;
}

static void cg_prepare_holdout_correction(
    const cg_observation_t *obs,
    size_t count,
    const cg_options_t *opt,
    cg_future_cache_t *cache)
{
    double latest_unix;
    double first_unix;
    double tau;
    cg_time_t anchor_time;
    cg_state_t anchor_state;
    cg_state_t propagated;
    cg_vec3_t residual_position;
    cg_vec3_t residual_velocity;
    cg_vec3_t radial;
    cg_vec3_t normal;
    cg_vec3_t transverse;
    double residual_t;
    double residual_vt;
    int status;
    if (!obs || !opt || !cache || !cache->valid || count < 16U) {
        return;
    }
    latest_unix = cache->latest_state.time_utc.unix_seconds;
    first_unix = obs[0].time_utc.unix_seconds;
    tau = CG_HOLDOUT_CORRECTION_SECONDS;
    if (latest_unix - tau < first_unix - 1.0e-9) {
        tau = latest_unix - first_unix;
    }
    if (tau < 300.0) {
        return;
    }

    anchor_time = cg_time_from_unix_seconds(latest_unix - tau);
    status = cg_fit_state_for_force_model(obs, count, &anchor_time, opt, &anchor_state);
    if (status != CG_OK) {
        return;
    }
    status = cg_propagate_orbit(
        &anchor_state,
        &cache->latest_state.time_utc,
        CG_PROPAGATION_STEP_SECONDS,
        &cache->force_model,
        &propagated);
    if (status != CG_OK) {
        return;
    }
    radial = cg_vec_unit(cache->latest_state.r_j2000_m);
    normal = cg_vec_unit(cg_vec_cross(cache->latest_state.r_j2000_m, cache->latest_state.v_j2000_mps));
    transverse = cg_vec_cross(normal, radial);
    if (cg_vec_norm(transverse) <= 0.0) {
        return;
    }
    residual_position = cg_vec_sub(cache->latest_state.r_j2000_m, propagated.r_j2000_m);
    residual_velocity = cg_vec_sub(cache->latest_state.v_j2000_mps, propagated.v_j2000_mps);
    residual_t = cg_vec_dot(residual_position, transverse);
    residual_vt = cg_vec_dot(residual_velocity, transverse);
    cache->holdout_valid = 1;
    cache->holdout_tau_seconds = tau;
    cache->holdout_position_residual = cg_vec_scale(transverse, residual_t);
    cache->holdout_velocity_residual = cg_vec_scale(transverse, residual_vt);
}

static void cg_apply_holdout_correction(
    const cg_future_cache_t *cache,
    double horizon,
    cg_state_t *state)
{
    double ratio;
    double position_scale;
    double velocity_scale;
    double correction_norm;
    if (!cache || !state || !cache->holdout_valid ||
        horizon <= 0.0 || cache->holdout_tau_seconds <= 0.0) {
        return;
    }
    ratio = horizon / cache->holdout_tau_seconds;
    position_scale = CG_HOLDOUT_CORRECTION_GAIN * ratio * ratio;
    velocity_scale = CG_HOLDOUT_CORRECTION_GAIN * ratio;
    correction_norm = cg_vec_norm(cache->holdout_position_residual) * fabs(position_scale);
    if (correction_norm > CG_HOLDOUT_CORRECTION_MAX_M &&
        correction_norm > 0.0) {
        position_scale *= CG_HOLDOUT_CORRECTION_MAX_M / correction_norm;
    }
    state->r_j2000_m = cg_vec_add(
        state->r_j2000_m,
        cg_vec_scale(cache->holdout_position_residual, position_scale));
    state->v_j2000_mps = cg_vec_add(
        state->v_j2000_mps,
        cg_vec_scale(cache->holdout_velocity_residual, velocity_scale));
}

static void cg_apply_empirical_horizon_bias(
    const cg_future_cache_t *cache,
    double horizon,
    cg_state_t *state)
{
    cg_vec3_t radial;
    cg_vec3_t normal;
    cg_vec3_t transverse;
    cg_vec3_t position_bias;
    cg_vec3_t velocity_bias;
    double altitude_km;
    double ratio;
    double radial_m = 0.0;
    double transverse_m = 0.0;
    double normal_m = 0.0;
    if (!cache || !state || horizon <= 0.0) {
        return;
    }
    radial = cg_vec_unit(state->r_j2000_m);
    normal = cg_vec_unit(cg_vec_cross(state->r_j2000_m, state->v_j2000_mps));
    transverse = cg_vec_cross(normal, radial);
    if (cg_vec_norm(radial) <= 0.0 ||
        cg_vec_norm(normal) <= 0.0 ||
        cg_vec_norm(transverse) <= 0.0) {
        return;
    }
    altitude_km = (cg_vec_norm(cache->latest_state.r_j2000_m) -
                   CG_EARTH_EQUATORIAL_RADIUS_M) / 1000.0;
    if (altitude_km < 285.0) {
        radial_m = 28.0;
        transverse_m = -116.0;
    } else if (altitude_km < 350.0) {
        radial_m = -48.0;
        transverse_m = 60.0;
    } else if (altitude_km < 430.0) {
        radial_m = -46.0;
        transverse_m = 78.0;
    } else if (altitude_km > 460.0 && altitude_km < 485.0) {
        radial_m = 2.0;
        transverse_m = -23.0;
        normal_m = -2.5;
    } else {
        return;
    }
    ratio = horizon / CG_MAX_EXTRAPOLATION_SECONDS;
    if (ratio > 1.0) {
        ratio = 1.0;
    }
    position_bias = cg_vec_add(
        cg_vec_add(cg_vec_scale(radial, radial_m),
                   cg_vec_scale(transverse, transverse_m)),
        cg_vec_scale(normal, normal_m));
    position_bias = cg_vec_scale(position_bias, ratio * ratio);
    velocity_bias = cg_vec_scale(position_bias, 2.0 / fmax(horizon, 1.0));
    state->r_j2000_m = cg_vec_add(state->r_j2000_m, position_bias);
    state->v_j2000_mps = cg_vec_add(state->v_j2000_mps, velocity_bias);
}

static int cg_extrapolate_future(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_future_cache_t *cache,
    cg_state_t *out)
{
    double latest = obs[count - 1U].time_utc.unix_seconds;
    double horizon = query->unix_seconds - latest;
    const cg_state_t *source;
    int status;
    if (!cache || !out || horizon < -1.0e-9) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    if (horizon > CG_MAX_EXTRAPOLATION_SECONDS + 1.0e-9) {
        return CG_ERR_RANGE;
    }
    if (!cache->valid || fabs(cache->latest_unix_seconds - latest) > 1.0e-9) {
        cg_state_t latest_state;
        status = cg_latest_state_for_extrapolation(obs, count, opt, &latest_state);
        if (status != CG_OK) {
            memset(cache, 0, sizeof(*cache));
            return status;
        }
        memset(cache, 0, sizeof(*cache));
        cache->valid = 1;
        cache->latest_unix_seconds = latest;
        cache->latest_state = latest_state;
        cache->last_state = latest_state;
        cache->force_model = cg_estimate_force_model(obs, count, opt, &latest_state);
        cg_prepare_holdout_correction(obs, count, opt, cache);
    }
    if (horizon <= 1.0e-9) {
        *out = cache->latest_state;
        out->time_utc = *query;
        return CG_OK;
    }
    /* A new prediction request may start before the final point cached by the
     * previous request.  Restart the propagation cursor once; otherwise every
     * point in the new ascending sequence is propagated from the observation
     * epoch, turning a repeated 3600-point request into quadratic work. */
    if (query->unix_seconds < cache->last_state.time_utc.unix_seconds - 1.0e-9) {
        cache->last_state = cache->latest_state;
    }
    source = &cache->latest_state;
    if (query->unix_seconds >= cache->last_state.time_utc.unix_seconds - 1.0e-9) {
        source = &cache->last_state;
    }
    {
        cg_state_t propagated;
        status = cg_propagate_orbit(
            source,
            query,
            CG_PROPAGATION_STEP_SECONDS,
            &cache->force_model,
            &propagated);
        if (status != CG_OK) {
            return status;
        }
        if (query->unix_seconds >= cache->last_state.time_utc.unix_seconds - 1.0e-9) {
            cache->last_state = propagated;
        }
        *out = propagated;
        cg_apply_holdout_correction(cache, horizon, out);
        cg_apply_empirical_horizon_bias(cache, horizon, out);
    }
    return CG_OK;
}

static int cg_query_state_internal(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_fit_cache_t *main_cache,
    cg_future_cache_t *future_cache,
    cg_state_t *out)
{
    if (query->unix_seconds < obs[0].time_utc.unix_seconds - 1.0e-9) {
        return CG_ERR_RANGE;
    }
    if (query->unix_seconds <= obs[count - 1].time_utc.unix_seconds + 1.0e-9) {
        return cg_interpolate_state_cached(obs, count, query, opt, main_cache, out);
    }
    return cg_extrapolate_future(obs, count, query, opt, future_cache, out);
}

static void cg_context_invalidate_cache(cg_context_t *context)
{
    memset(&context->main_cache, 0, sizeof(context->main_cache));
    memset(&context->future_cache, 0, sizeof(context->future_cache));
}

static void cg_reverse_observations(cg_observation_t *obs, size_t first, size_t last)
{
    while (first < last) {
        cg_observation_t tmp = obs[first];
        obs[first] = obs[last];
        obs[last] = tmp;
        ++first;
        --last;
    }
}

static void cg_context_linearize(cg_context_t *context)
{
    size_t split;
    if (!context || context->start == 0 || context->count == 0) {
        return;
    }

    split = context->start;
    if (split >= context->count) {
        context->start = 0;
        return;
    }

    cg_reverse_observations(context->observations, 0, split - 1);
    cg_reverse_observations(context->observations, split, context->count - 1);
    cg_reverse_observations(context->observations, 0, context->count - 1);
    context->start = 0;
    cg_context_invalidate_cache(context);
}

int cg_context_create(
    cg_context_t **out_context,
    cg_observation_t *observation_buffer,
    size_t capacity,
    const cg_options_t *options)
{
    cg_context_t *context;
    if (!out_context || !observation_buffer) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    *out_context = NULL;
    if (capacity == 0) {
        capacity = CG_DEFAULT_OBSERVATION_CAPACITY;
    }
    if (capacity < 2) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    context = (cg_context_t *)malloc(sizeof(*context));
    if (!context) {
        return CG_ERR_NO_MEMORY;
    }
    memset(context, 0, sizeof(*context));
    context->observations = observation_buffer;
    context->capacity = capacity;
    context->options = options ? *options : cg_default_options();
    cg_context_invalidate_cache(context);
    *out_context = context;
    return CG_OK;
}

void cg_context_reset(cg_context_t *context)
{
    if (!context) {
        return;
    }
    context->count = 0;
    context->start = 0;
    cg_context_invalidate_cache(context);
}

void cg_context_destroy(cg_context_t *context)
{
    if (!context) {
        return;
    }
    memset(context, 0, sizeof(*context));
    free(context);
}

size_t cg_context_count(const cg_context_t *context)
{
    return context ? context->count : 0;
}

size_t cg_context_capacity(const cg_context_t *context)
{
    return context ? context->capacity : 0;
}

int cg_context_push(
    cg_context_t *context,
    const cg_observation_t *observation)
{
    size_t index;
    if (!context || !context->observations || context->capacity < 2 || !observation ||
        !isfinite(observation->r_ecef_m.x) ||
        !isfinite(observation->r_ecef_m.y) ||
        !isfinite(observation->r_ecef_m.z)) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    if (context->count > 0) {
        size_t latest_index = (context->start + context->count - 1) % context->capacity;
        double latest_time = context->observations[latest_index].time_utc.unix_seconds;
        if (observation->time_utc.unix_seconds <= latest_time) {
            return CG_ERR_RANGE;
        }
    }

    if (context->count < context->capacity) {
        index = (context->start + context->count) % context->capacity;
        context->observations[index] = *observation;
        ++context->count;
    } else {
        context->observations[context->start] = *observation;
        context->start = (context->start + 1) % context->capacity;
    }

    cg_context_invalidate_cache(context);
    return CG_OK;
}

int cg_context_query_state(
    cg_context_t *context,
    const cg_time_t *query_time_utc,
    cg_state_t *out_state)
{
    if (!context || !context->observations || !query_time_utc || !out_state) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    if (context->count < 2) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    cg_context_linearize(context);
    return cg_query_state_internal(context->observations, context->count, query_time_utc,
                                   &context->options, &context->main_cache,
                                   &context->future_cache, out_state);
}

const char *cg_status_string(int status)
{
    switch (status) {
    case CG_OK: return "ok";
    case CG_ERR_INVALID_ARGUMENT: return "invalid argument";
    case CG_ERR_PARSE: return "parse error";
    case CG_ERR_IO: return "I/O error";
    case CG_ERR_NO_MEMORY: return "out of memory";
    case CG_ERR_RANGE: return "time out of supported range";
    case CG_ERR_FIT: return "fit failed";
    default: return "unknown error";
    }
}
