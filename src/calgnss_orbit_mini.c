#include "calgnss_orbit_mini.h"

#include <float.h>
#include <math.h>
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

struct cg_context_t {
    cg_observation_t *observations;
    size_t capacity;
    size_t count;
    size_t start;
    cg_options_t options;
    cg_fit_cache_t main_cache;
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
    size_t count = sizeof(k_cg_eop_2026) / sizeof(k_cg_eop_2026[0]);
    size_t i;
    cg_eop_value_t empty;
    empty.dut1_seconds = 0.0;
    empty.x_pole_rad = 0.0;
    empty.y_pole_rad = 0.0;
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

static int cg_query_state_internal(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_fit_cache_t *main_cache,
    cg_state_t *out)
{
    if (query->unix_seconds < obs[0].time_utc.unix_seconds - 1.0e-9) {
        return CG_ERR_RANGE;
    }
    if (query->unix_seconds <= obs[count - 1].time_utc.unix_seconds + 1.0e-9) {
        return cg_interpolate_state_cached(obs, count, query, opt, main_cache, out);
    }
    return CG_ERR_RANGE;
}

static void cg_context_invalidate_cache(cg_context_t *context)
{
    memset(&context->main_cache, 0, sizeof(context->main_cache));
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
                                   &context->options, &context->main_cache, out_state);
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
