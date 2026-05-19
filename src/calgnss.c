#include "calgnss.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifndef CG_PI
#define CG_PI 3.141592653589793238462643383279502884
#endif

#define CG_D2PI (2.0 * CG_PI)
#define CG_DEG2RAD (CG_PI / 180.0)
#define CG_ARCSEC2RAD (CG_PI / (180.0 * 3600.0))
#define CG_MAS2RAD (CG_ARCSEC2RAD / 1000.0)
#define CG_JD_J2000 2451545.0
#define CG_JULIAN_CENTURY 36525.0
#define CG_TT_MINUS_UTC_SECONDS 69.184

#define CG_MU_EARTH 3.986004418e14
#define CG_OMEGA_EARTH 7.2921150e-5
#define CG_J2_EARTH 1.08262668e-3
#define CG_WGS84_A 6378137.0
#define CG_WGS84_F (1.0 / 298.257223563)
#define CG_MAX_DEGREE 16
#define CG_DIRECT_VELOCITY_WINDOW 9
#define CG_DIRECT_VELOCITY_DEGREE 8

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

static double cg_vec_dot(cg_vec3_t a, cg_vec3_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

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

static cg_vec3_t cg_vec_scale(cg_vec3_t a, double s)
{
    cg_vec3_t r;
    r.x = a.x * s;
    r.y = a.y * s;
    r.z = a.z * s;
    return r;
}

static double cg_vec_norm(cg_vec3_t a)
{
    return sqrt(cg_vec_dot(a, a));
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

static cg_time_t cg_time_from_calendar(int y, int m, int d, int hour, int minute, double second)
{
    cg_time_t t;
    int whole_second = (int)floor(second);
    double frac = second - (double)whole_second;
    int64_t days = cg_days_from_civil(y, (unsigned)m, (unsigned)d);
    t.year = y;
    t.month = m;
    t.day = d;
    t.hour = hour;
    t.minute = minute;
    t.second = second;
    t.jd_utc = cg_jd_from_calendar(y, m, d, hour, minute, second);
    t.unix_seconds = (double)days * 86400.0 + (double)hour * 3600.0 +
                     (double)minute * 60.0 + (double)whole_second + frac;
    return t;
}

static cg_time_t cg_time_add_seconds(const cg_time_t *time_utc, double seconds)
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
    cg_civil_from_days(days, &y, &m, &d);
    hour = (int)(sod / 3600);
    minute = (int)((sod % 3600) / 60);
    return cg_time_from_calendar(y, (int)m, (int)d, hour, minute,
                                 (double)(sod % 60) + frac);
}

static cg_vec3_t cg_geodetic_to_ecef(double lat_deg, double lon_deg, double h_m)
{
    const double e2 = CG_WGS84_F * (2.0 - CG_WGS84_F);
    double lat = lat_deg * CG_DEG2RAD;
    double lon = lon_deg * CG_DEG2RAD;
    double sin_lat = sin(lat);
    double cos_lat = cos(lat);
    double sin_lon = sin(lon);
    double cos_lon = cos(lon);
    double n = CG_WGS84_A / sqrt(1.0 - e2 * sin_lat * sin_lat);
    cg_vec3_t r;
    r.x = (n + h_m) * cos_lat * cos_lon;
    r.y = (n + h_m) * cos_lat * sin_lon;
    r.z = (n * (1.0 - e2) + h_m) * sin_lat;
    return r;
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
    double jd_tt = jd_utc + CG_TT_MINUS_UTC_SECONDS / 86400.0;
    double rbpn[3][3];
    double rbpn_t[3][3];
    double rz[3][3];
    double dpsi;
    double epsa;
    double gst;
    int i;
    int j;
    cg_pnm00b(jd_tt, rbpn, &dpsi, &epsa);
    gst = cg_anp(cg_gmst00(jd_utc, jd_tt) + dpsi * cos(epsa));
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            rbpn_t[i][j] = rbpn[j][i];
        }
    }
    cg_rz_vector_matrix(gst, rz);
    cg_mat_mul(rbpn_t, rz, out);
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
    cg_vec3_t r_ecef;
    cg_vec3_t r_j2000;
    cg_vec3_t ignored_v_j2000;
    zero.x = zero.y = zero.z = 0.0;
    r_ecef = cg_geodetic_to_ecef(
        observation->lat_deg, observation->lon_deg, observation->alt_m);
    cg_ecef_to_j2000(
        &observation->time_utc,
        r_ecef,
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
    opt.fit_window_minutes = 60.0;
    opt.max_extrapolation_seconds = 0.0;
    opt.propagation_step_seconds = 10.0;
    opt.enable_orbit_phase_correction = 1;
    return opt;
}

static void cg_select_window(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    size_t *first,
    size_t *last)
{
    double half_window = opt->fit_window_minutes * 30.0;
    double q = query->unix_seconds;
    size_t i;
    size_t f = 0;
    size_t l = count - 1;
    if (opt->fit_window_minutes <= 0.0) {
        *first = 0;
        *last = count - 1;
        return;
    }
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

static int cg_estimate_velocity_at_index(
    const cg_observation_t *obs,
    size_t count,
    size_t index,
    cg_vec3_t *out_velocity)
{
    size_t first;
    size_t last;
    size_t nobs;
    size_t i;
    int degree;
    int ncoef;
    int j;
    int k;
    int coord;
    double center_seconds;
    double scale_seconds = 0.0;
    double ata[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 1];
    double rhs[3][CG_MAX_DEGREE + 1];

    if (!obs || !out_velocity || index >= count || count < 2) {
        return CG_ERR_INVALID_ARGUMENT;
    }

    first = index;
    if (first > CG_DIRECT_VELOCITY_WINDOW / 2) {
        first -= CG_DIRECT_VELOCITY_WINDOW / 2;
    } else {
        first = 0;
    }
    last = first + CG_DIRECT_VELOCITY_WINDOW - 1;
    if (last >= count) {
        last = count - 1;
        first = last + 1 > CG_DIRECT_VELOCITY_WINDOW ? last + 1 - CG_DIRECT_VELOCITY_WINDOW : 0;
    }
    nobs = last - first + 1;
    degree = (int)nobs - 1;
    if (degree > CG_DIRECT_VELOCITY_DEGREE) {
        degree = CG_DIRECT_VELOCITY_DEGREE;
    }
    if (degree < 1) {
        return CG_ERR_FIT;
    }
    ncoef = degree + 1;
    center_seconds = obs[index].time_utc.unix_seconds;
    for (i = first; i <= last; ++i) {
        double dt = fabs(obs[i].time_utc.unix_seconds - center_seconds);
        if (dt > scale_seconds) {
            scale_seconds = dt;
        }
    }
    if (scale_seconds <= 0.0) {
        return CG_ERR_FIT;
    }

    memset(ata, 0, sizeof(ata));
    memset(rhs, 0, sizeof(rhs));
    for (i = first; i <= last; ++i) {
        double powers[CG_MAX_DEGREE + 1];
        double tau = (obs[i].time_utc.unix_seconds - center_seconds) / scale_seconds;
        cg_vec3_t r_j2000 = cg_observation_to_j2000_position(&obs[i]);
        powers[0] = 1.0;
        for (j = 1; j < ncoef; ++j) {
            powers[j] = powers[j - 1] * tau;
        }
        for (j = 0; j < ncoef; ++j) {
            for (k = 0; k < ncoef; ++k) {
                ata[j][k] += powers[j] * powers[k];
            }
            rhs[0][j] += powers[j] * r_j2000.x;
            rhs[1][j] += powers[j] * r_j2000.y;
            rhs[2][j] += powers[j] * r_j2000.z;
        }
    }

    out_velocity->x = 0.0;
    out_velocity->y = 0.0;
    out_velocity->z = 0.0;
    for (coord = 0; coord < 3; ++coord) {
        double aug[CG_MAX_DEGREE + 1][CG_MAX_DEGREE + 2];
        double coeff[CG_MAX_DEGREE + 1];
        int status;
        memset(aug, 0, sizeof(aug));
        memset(coeff, 0, sizeof(coeff));
        for (j = 0; j < ncoef; ++j) {
            for (k = 0; k < ncoef; ++k) {
                aug[j][k] = ata[j][k];
            }
            aug[j][ncoef] = rhs[coord][j];
        }
        status = cg_solve_linear(ncoef, aug, coeff);
        if (status != CG_OK) {
            return status;
        }
        if (coord == 0) {
            out_velocity->x = coeff[1] / scale_seconds;
        } else if (coord == 1) {
            out_velocity->y = coeff[1] / scale_seconds;
        } else {
            out_velocity->z = coeff[1] / scale_seconds;
        }
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
    cg_select_window(obs, count, query, opt, &first, &last);
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

static cg_vec3_t cg_j2_acceleration(cg_vec3_t r)
{
    double radius = cg_vec_norm(r);
    double r2 = radius * radius;
    double z2 = r.z * r.z;
    double scale = 1.5 * CG_J2_EARTH * CG_MU_EARTH * CG_WGS84_A * CG_WGS84_A / pow(radius, 5.0);
    double z_factor = 5.0 * z2 / r2;
    cg_vec3_t a;
    double two = -CG_MU_EARTH / (radius * radius * radius);
    a.x = two * r.x + scale * r.x * (z_factor - 1.0);
    a.y = two * r.y + scale * r.y * (z_factor - 1.0);
    a.z = two * r.z + scale * r.z * (z_factor - 3.0);
    return a;
}

static void cg_state_derivative(const double y[6], double dy[6])
{
    cg_vec3_t r;
    cg_vec3_t a;
    r.x = y[0]; r.y = y[1]; r.z = y[2];
    a = cg_j2_acceleration(r);
    dy[0] = y[3]; dy[1] = y[4]; dy[2] = y[5];
    dy[3] = a.x; dy[4] = a.y; dy[5] = a.z;
}

static int cg_propagate_j2(cg_vec3_t r0, cg_vec3_t v0, double dt, double step_seconds, cg_vec3_t *r, cg_vec3_t *v)
{
    double y[6];
    double remaining;
    double direction;
    int i;
    if (step_seconds <= 0.0) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    y[0] = r0.x; y[1] = r0.y; y[2] = r0.z;
    y[3] = v0.x; y[4] = v0.y; y[5] = v0.z;
    remaining = fabs(dt);
    direction = dt >= 0.0 ? 1.0 : -1.0;
    while (remaining > 1.0e-12) {
        double h = direction * fmin(step_seconds, remaining);
        double k1[6], k2[6], k3[6], k4[6], tmp[6];
        cg_state_derivative(y, k1);
        for (i = 0; i < 6; ++i) tmp[i] = y[i] + 0.5 * h * k1[i];
        cg_state_derivative(tmp, k2);
        for (i = 0; i < 6; ++i) tmp[i] = y[i] + 0.5 * h * k2[i];
        cg_state_derivative(tmp, k3);
        for (i = 0; i < 6; ++i) tmp[i] = y[i] + h * k3[i];
        cg_state_derivative(tmp, k4);
        for (i = 0; i < 6; ++i) {
            y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        remaining -= fabs(h);
    }
    r->x = y[0]; r->y = y[1]; r->z = y[2];
    v->x = y[3]; v->y = y[4]; v->z = y[5];
    return CG_OK;
}

static int cg_direct_state_at_index(
    const cg_observation_t *obs,
    size_t count,
    size_t index,
    const cg_options_t *opt,
    cg_fit_cache_t *cache,
    cg_state_t *out)
{
    int status;
    if (index >= count) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    out->time_utc = obs[index].time_utc;
    out->r_j2000_m = cg_observation_to_j2000_position(&obs[index]);
    status = cg_estimate_velocity_at_index(obs, count, index, &out->v_j2000_mps);
    if (status == CG_OK) {
        return CG_OK;
    }
    status = cg_interpolate_state_cached(obs, count, &obs[index].time_utc, opt, cache, out);
    if (status == CG_OK) {
        out->r_j2000_m = cg_observation_to_j2000_position(&obs[index]);
    }
    return status;
}

static size_t cg_find_nearest_observation(const cg_observation_t *obs, size_t count, double unix_seconds)
{
    size_t lo = 0;
    size_t hi = count;
    if (unix_seconds <= obs[0].time_utc.unix_seconds) {
        return 0;
    }
    if (unix_seconds >= obs[count - 1].time_utc.unix_seconds) {
        return count - 1;
    }
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (obs[mid].time_utc.unix_seconds < unix_seconds) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo == 0) {
        return 0;
    }
    if (fabs(obs[lo].time_utc.unix_seconds - unix_seconds) <
        fabs(obs[lo - 1].time_utc.unix_seconds - unix_seconds)) {
        return lo;
    }
    return lo - 1;
}

static double cg_orbital_period_seconds(cg_vec3_t r, cg_vec3_t v)
{
    double rn = cg_vec_norm(r);
    double v2 = cg_vec_dot(v, v);
    double denom = 2.0 / rn - v2 / CG_MU_EARTH;
    double a;
    if (denom <= 0.0) {
        return 0.0;
    }
    a = 1.0 / denom;
    if (a <= 0.0) {
        return 0.0;
    }
    return CG_D2PI * sqrt((a * a * a) / CG_MU_EARTH);
}

static int cg_extrapolate_future(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_fit_cache_t *cache,
    cg_state_t *out)
{
    double horizon = query->unix_seconds - obs[count - 1].time_utc.unix_seconds;
    cg_state_t end_state;
    cg_vec3_t pred_r;
    cg_vec3_t pred_v;
    int status;
    if (horizon < -1.0e-9) {
        return CG_ERR_RANGE;
    }
    if (opt->max_extrapolation_seconds > 0.0 &&
        horizon > opt->max_extrapolation_seconds + 1.0e-9) {
        return CG_ERR_RANGE;
    }
    status = cg_direct_state_at_index(obs, count, count - 1, opt, cache, &end_state);
    if (status != CG_OK) {
        return status;
    }
    status = cg_propagate_j2(end_state.r_j2000_m, end_state.v_j2000_mps, horizon,
                             opt->propagation_step_seconds, &pred_r, &pred_v);
    if (status != CG_OK) {
        return status;
    }
    if (opt->enable_orbit_phase_correction && horizon > 0.0) {
        double period = cg_orbital_period_seconds(end_state.r_j2000_m, end_state.v_j2000_mps);
        if (period > horizon + 60.0) {
            double anchor_unix = end_state.time_utc.unix_seconds - period;
            size_t anchor_index = cg_find_nearest_observation(obs, count, anchor_unix);
            cg_state_t anchor_state;
            cg_state_t target_state;
            cg_time_t target_time;
            cg_vec3_t past_pred_r;
            cg_vec3_t past_pred_v;
            double actual_horizon;
            memset(&anchor_state, 0, sizeof(anchor_state));
            memset(&target_state, 0, sizeof(target_state));
            memset(&past_pred_r, 0, sizeof(past_pred_r));
            memset(&past_pred_v, 0, sizeof(past_pred_v));
            if (anchor_index < count - 1 &&
                fabs(obs[anchor_index].time_utc.unix_seconds - anchor_unix) <= 180.0) {
                target_time = cg_time_add_seconds(&obs[anchor_index].time_utc, horizon);
                if (target_time.unix_seconds <= obs[count - 1].time_utc.unix_seconds + 1.0e-9 &&
                    target_time.unix_seconds >= obs[0].time_utc.unix_seconds - 1.0e-9) {
                    status = cg_direct_state_at_index(obs, count, anchor_index, opt, cache, &anchor_state);
                    if (status == CG_OK) {
                        status = cg_interpolate_state_cached(obs, count, &target_time, opt, cache, &target_state);
                    }
                    if (status == CG_OK) {
                        actual_horizon = target_time.unix_seconds - anchor_state.time_utc.unix_seconds;
                        status = cg_propagate_j2(anchor_state.r_j2000_m, anchor_state.v_j2000_mps,
                                                 actual_horizon, opt->propagation_step_seconds,
                                                 &past_pred_r, &past_pred_v);
                    }
                    if (status == CG_OK) {
                        cg_vec3_t correction = cg_vec_sub(target_state.r_j2000_m, past_pred_r);
                        pred_r = cg_vec_add(pred_r, correction);
                    }
                }
            }
        }
    }
    out->time_utc = *query;
    out->r_j2000_m = pred_r;
    out->v_j2000_mps = pred_v;
    return CG_OK;
}

static int cg_query_state_internal(
    const cg_observation_t *obs,
    size_t count,
    const cg_time_t *query,
    const cg_options_t *opt,
    cg_fit_cache_t *main_cache,
    cg_fit_cache_t *aux_cache,
    cg_state_t *out)
{
    if (query->unix_seconds < obs[0].time_utc.unix_seconds - 1.0e-9) {
        return CG_ERR_RANGE;
    }
    if (query->unix_seconds <= obs[count - 1].time_utc.unix_seconds + 1.0e-9) {
        return cg_interpolate_state_cached(obs, count, query, opt, main_cache, out);
    }
    return cg_extrapolate_future(obs, count, query, opt, aux_cache, out);
}

int cg_query_state(
    const cg_observation_t *observations,
    size_t count,
    const cg_time_t *query_time_utc,
    const cg_options_t *options,
    cg_state_t *out_state)
{
    cg_options_t default_options;
    cg_fit_cache_t main_cache;
    cg_fit_cache_t aux_cache;
    if (!observations || count < 2 || !query_time_utc || !out_state) {
        return CG_ERR_INVALID_ARGUMENT;
    }
    default_options = options ? *options : cg_default_options();
    memset(&main_cache, 0, sizeof(main_cache));
    memset(&aux_cache, 0, sizeof(aux_cache));
    return cg_query_state_internal(observations, count, query_time_utc, &default_options,
                                   &main_cache, &aux_cache, out_state);
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
