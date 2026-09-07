/* orbit_determination.c -- standalone C99 orbit determination implementation.
 * Input: RTCM3 / NMEA. Output: latest observed J2000 position and velocity.
 * All implementation dependencies and EOP data are embedded below.
 * RTKLIB copyright/license: see embedded notice and third_party/rtklib/LICENSE.txt.
 */
#ifndef _DEFAULT_SOURCE
#define _DEFAULT_SOURCE
#endif
#include "orbit_determination.h"
#include <errno.h>
#include <limits.h>
/* Extracted from PredictOrbit: coordinate conversion and fitting only. */
typedef struct od_fit_time_t {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    double jd_utc;
    double unix_seconds;
} od_fit_time_t;

typedef struct od_fit_vec3_t {
    double x;
    double y;
    double z;
} od_fit_vec3_t;

typedef struct od_fit_state_t {
    od_fit_time_t time_utc;
    od_fit_vec3_t r_j2000_m;
    od_fit_vec3_t v_j2000_mps;
} od_fit_state_t;

typedef struct od_fit_observation_t {
    od_fit_time_t time_utc;
    od_fit_vec3_t r_ecef_m;
} od_fit_observation_t;

typedef struct od_fit_options_t {
    int degree;
} od_fit_options_t;

#define OD_FIT_MAX_DEGREE 16
#define OD_FIT_DEFAULT_OBSERVATION_CAPACITY 10

typedef struct od_fit_context_t od_fit_context_t;

typedef enum od_fit_status_t {
    OD_FIT_OK = 0,
    OD_FIT_ERR_INVALID_ARGUMENT = -1,
    OD_FIT_ERR_PARSE = -2,
    OD_FIT_ERR_IO = -3,
    OD_FIT_ERR_NO_MEMORY = -4,
    OD_FIT_ERR_RANGE = -5,
    OD_FIT_ERR_FIT = -6
} od_fit_status_t;

static od_fit_options_t od_fit_default_options(void);


static int od_fit_context_create(
    od_fit_context_t **out_context,
    od_fit_observation_t *observation_buffer,
    size_t capacity,
    const od_fit_options_t *options);

static void od_fit_context_reset(od_fit_context_t *context);

static void od_fit_context_destroy(od_fit_context_t *context);

static size_t od_fit_context_count(const od_fit_context_t *context);


static int od_fit_context_push(
    od_fit_context_t *context,
    const od_fit_observation_t *observation);

static int od_fit_context_query_state(
    od_fit_context_t *context,
    const od_fit_time_t *query_time_utc,
    od_fit_state_t *out_state);




#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef OD_FIT_PI
#define OD_FIT_PI 3.141592653589793238462643383279502884
#endif

#define OD_FIT_D2PI (2.0 * OD_FIT_PI)
#define OD_FIT_ARCSEC2RAD (OD_FIT_PI / (180.0 * 3600.0))
#define OD_FIT_MAS2RAD (OD_FIT_ARCSEC2RAD / 1000.0)
#define OD_FIT_JD_J2000 2451545.0
#define OD_FIT_JULIAN_CENTURY 36525.0
#define OD_FIT_TT_MINUS_UTC_SECONDS 69.184
#define OD_FIT_SECONDS_PER_DAY 86400.0

#define OD_FIT_OMEGA_EARTH 7.2921150e-5

#ifndef OD_FIT_USE_EOP_TABLE
#define OD_FIT_USE_EOP_TABLE 1
#endif

typedef struct od_fit_fit_t {
    int degree;
    size_t first_index;
    size_t last_index;
    double epoch_seconds;
    double scale_seconds;
    double coeff[3][OD_FIT_MAX_DEGREE + 1];
} od_fit_fit_t;

typedef struct od_fit_fit_cache_t {
    int valid;
    int degree;
    size_t first_index;
    size_t last_index;
    od_fit_fit_t fit;
} od_fit_fit_cache_t;

struct od_fit_context_t {
    od_fit_observation_t *observations;
    size_t capacity;
    size_t count;
    size_t start;
    od_fit_options_t options;
    od_fit_fit_cache_t main_cache;
};

typedef struct od_fit_eop_record_t {
    int mjd;
    double dut1_seconds;
    double x_pole_arcsec;
    double y_pole_arcsec;
} od_fit_eop_record_t;

typedef struct od_fit_eop_value_t {
    double dut1_seconds;
    double x_pole_rad;
    double y_pole_rad;
} od_fit_eop_value_t;

// Generated from IERS finals2000A.data.csv. Units: MJD, UT1-UTC seconds, pole arcseconds.
static const od_fit_eop_record_t od_eop_2026[] = {
    {61041, 0.0740684, 0.110516, 0.331199},
    {61042, 0.0741649, 0.109622, 0.332518},
    {61043, 0.0743605, 0.108264, 0.333542},
    {61044, 0.0744805, 0.106883, 0.334545},
    {61045, 0.0743581, 0.105523, 0.335519},
    {61046, 0.0740241, 0.104038, 0.336285},
    {61047, 0.0734710, 0.102511, 0.336912},
    {61048, 0.0728581, 0.101216, 0.337391},
    {61049, 0.0722579, 0.100610, 0.337813},
    {61050, 0.0717709, 0.100539, 0.338547},
    {61051, 0.0714303, 0.100522, 0.339366},
    {61052, 0.0713231, 0.100288, 0.340025},
    {61053, 0.0714171, 0.099850, 0.340548},
    {61054, 0.0716604, 0.099338, 0.340914},
    {61055, 0.0721092, 0.098524, 0.341279},
    {61056, 0.0726595, 0.097299, 0.341630},
    {61057, 0.0732641, 0.096058, 0.341985},
    {61058, 0.0738300, 0.094802, 0.342574},
    {61059, 0.0742076, 0.093417, 0.343374},
    {61060, 0.0743258, 0.092272, 0.344529},
    {61061, 0.0741780, 0.091599, 0.345989},
    {61062, 0.0737311, 0.091246, 0.347616},
    {61063, 0.0730554, 0.091170, 0.349213},
    {61064, 0.0722133, 0.091562, 0.350709},
    {61065, 0.0714661, 0.091903, 0.352312},
    {61066, 0.0708885, 0.091880, 0.353652},
    {61067, 0.0705113, 0.091936, 0.354778},
    {61068, 0.0703699, 0.092033, 0.356060},
    {61069, 0.0703858, 0.091846, 0.357430},
    {61070, 0.0705721, 0.091660, 0.358790},
    {61071, 0.0706679, 0.091704, 0.360202},
    {61072, 0.0706631, 0.091850, 0.361473},
    {61073, 0.0703823, 0.091937, 0.362494},
    {61074, 0.0699105, 0.091931, 0.363434},
    {61075, 0.0692672, 0.092070, 0.364388},
    {61076, 0.0685727, 0.092337, 0.365664},
    {61077, 0.0679629, 0.092733, 0.367129},
    {61078, 0.0674779, 0.093226, 0.368366},
    {61079, 0.0671872, 0.093942, 0.369169},
    {61080, 0.0671420, 0.095303, 0.369998},
    {61081, 0.0672916, 0.096925, 0.371304},
    {61082, 0.0676495, 0.098453, 0.372814},
    {61083, 0.0681080, 0.099667, 0.374243},
    {61084, 0.0686723, 0.100557, 0.375484},
    {61085, 0.0691987, 0.101195, 0.376801},
    {61086, 0.0696466, 0.101493, 0.377810},
    {61087, 0.0699048, 0.101867, 0.378472},
    {61088, 0.0699031, 0.102533, 0.379508},
    {61089, 0.0695678, 0.103140, 0.380771},
    {61090, 0.0689784, 0.103865, 0.381623},
    {61091, 0.0682981, 0.104894, 0.382341},
    {61092, 0.0675845, 0.105913, 0.383409},
    {61093, 0.0670129, 0.106710, 0.384701},
    {61094, 0.0665920, 0.107286, 0.385853},
    {61095, 0.0663815, 0.107452, 0.386591},
    {61096, 0.0664264, 0.106894, 0.387212},
    {61097, 0.0667063, 0.106030, 0.387719},
    {61098, 0.0670399, 0.105602, 0.387937},
    {61099, 0.0672241, 0.105431, 0.388379},
    {61100, 0.0672047, 0.105057, 0.389086},
    {61101, 0.0669114, 0.104976, 0.389650},
    {61102, 0.0663791, 0.105346, 0.390183},
    {61103, 0.0656687, 0.105546, 0.391067},
    {61104, 0.0648521, 0.105215, 0.391940},
    {61105, 0.0641099, 0.104895, 0.392353},
    {61106, 0.0634991, 0.104774, 0.392696},
    {61107, 0.0631049, 0.104062, 0.393085},
    {61108, 0.0629240, 0.103189, 0.393389},
    {61109, 0.0627736, 0.102977, 0.393847},
    {61110, 0.0626790, 0.103132, 0.394325},
    {61111, 0.0626456, 0.103042, 0.394882},
    {61112, 0.0626599, 0.102823, 0.395412},
    {61113, 0.0625265, 0.102821, 0.396049},
    {61114, 0.0621771, 0.103234, 0.396822},
    {61115, 0.0616040, 0.104055, 0.397659},
    {61116, 0.0608242, 0.104708, 0.398667},
    {61117, 0.0598380, 0.105088, 0.399480},
    {61118, 0.0586910, 0.105648, 0.400134},
    {61119, 0.0575076, 0.106148, 0.400851},
    {61120, 0.0564627, 0.106485, 0.401570},
    {61121, 0.0556377, 0.107129, 0.402231},
    {61122, 0.0551609, 0.108161, 0.402918},
    {61123, 0.0549353, 0.109421, 0.403452},
    {61124, 0.0549038, 0.110782, 0.404019},
    {61125, 0.0548405, 0.112273, 0.404992},
    {61126, 0.0546925, 0.113962, 0.406177},
    {61127, 0.0543313, 0.116054, 0.407112},
    {61128, 0.0536891, 0.118652, 0.407966},
    {61129, 0.0528435, 0.120973, 0.408777},
    {61130, 0.0519275, 0.122779, 0.409108},
    {61131, 0.0509677, 0.124635, 0.409398},
    {61132, 0.0500286, 0.126386, 0.410186},
    {61133, 0.0491947, 0.127866, 0.411255},
    {61134, 0.0485661, 0.129558, 0.412047},
    {61135, 0.0481416, 0.131587, 0.412647},
    {61136, 0.0479219, 0.133572, 0.413263},
    {61137, 0.0478038, 0.135323, 0.413558},
    {61138, 0.0477549, 0.136926, 0.413586},
    {61139, 0.0477038, 0.138461, 0.413645},
    {61140, 0.0475817, 0.139835, 0.413663},
    {61141, 0.0472943, 0.141187, 0.413821},
    {61142, 0.0468328, 0.142367, 0.414201},
    {61143, 0.0461257, 0.143283, 0.414450},
    {61144, 0.0452243, 0.144253, 0.414598},
    {61145, 0.0440894, 0.145188, 0.414921},
    {61146, 0.0428445, 0.145846, 0.415193},
    {61147, 0.0415941, 0.146675, 0.415279},
    {61148, 0.0405425, 0.147586, 0.415371},
    {61149, 0.0397739, 0.148386, 0.415673},
    {61150, 0.0392791, 0.149233, 0.415991},
    {61151, 0.0389979, 0.150046, 0.416042},
    {61152, 0.0388379, 0.150904, 0.416307},
    {61153, 0.0386405, 0.151701, 0.416752},
    {61154, 0.0383147, 0.152543, 0.417189},
    {61155, 0.0377862, 0.153577, 0.417834},
    {61156, 0.0370733, 0.154571, 0.418566},
    {61157, 0.0361992, 0.155213, 0.419187},
    {61158, 0.0352416, 0.155436, 0.419390},
    {61159, 0.0343503, 0.155450, 0.418989},
    {61160, 0.0335543, 0.155509, 0.418133},
    {61161, 0.0329226, 0.155975, 0.417425},
    {61162, 0.0324955, 0.156819, 0.417231},
    {61163, 0.0322343, 0.158032, 0.416990},
    {61164, 0.0320650, 0.159687, 0.416403},
    {61165, 0.0320886, 0.161237, 0.416052},
    {61166, 0.0321136, 0.162310, 0.416121},
    {61167, 0.0321700, 0.163219, 0.416106},
    {61168, 0.0321971, 0.164193, 0.415724},
    {61169, 0.0320546, 0.164941, 0.415182},
    {61170, 0.0316934, 0.165786, 0.414395},
    {61171, 0.0309750, 0.167010, 0.413421},
    {61172, 0.0301511, 0.168109, 0.412477},
    {61173, 0.0291758, 0.169055, 0.411800},
    {61174, 0.0281273, 0.169846, 0.411637},
    {61175, 0.0271499, 0.170536, 0.411771},
    {61176, 0.0263823, 0.171246, 0.412117},
    {61177, 0.0259113, 0.171931, 0.412524},
    {61178, 0.0256978, 0.172687, 0.412907},
    {61179, 0.0256285, 0.173472, 0.413152},
    {61180, 0.0255213, 0.174174, 0.413047},
    {61181, 0.0252260, 0.174715, 0.412691},
    {61182, 0.0247054, 0.175350, 0.412217},
    {61183, 0.0239732, 0.175960, 0.411676},
    {61184, 0.0231082, 0.176659, 0.411165},
    {61185, 0.0222088, 0.177418, 0.410700},
    {61186, 0.0213516, 0.178226, 0.410333},
    {61187, 0.0205896, 0.179047, 0.410017},
    {61188, 0.0199742, 0.179848, 0.409721},
    {61189, 0.0195548, 0.180604, 0.409428},
    {61190, 0.0193552, 0.181346, 0.409119},
    {61191, 0.0193777, 0.182084, 0.408786},
    {61192, 0.0195620, 0.182815, 0.408407},
    {61193, 0.0198290, 0.183538, 0.407995},
    {61194, 0.0201241, 0.184255, 0.407553},
    {61195, 0.0203976, 0.184975, 0.407080},
    {61196, 0.0205800, 0.185697, 0.406587},
    {61197, 0.0206192, 0.186422, 0.406080},
    {61198, 0.0204924, 0.187154, 0.405568},
    {61199, 0.0201867, 0.187889, 0.405050},
    {61200, 0.0197540, 0.188621, 0.404526},
    {61201, 0.0192660, 0.189345, 0.403997},
    {61202, 0.0188353, 0.190059, 0.403457},
    {61203, 0.0185925, 0.190763, 0.402903},
    {61204, 0.0186340, 0.191457, 0.402333},
    {61205, 0.0189753, 0.192142, 0.401746},
    {61206, 0.0195421, 0.192820, 0.401143},
    {61207, 0.0201950, 0.193491, 0.400524},
    {61208, 0.0207768, 0.194154, 0.399891},
    {61209, 0.0211740, 0.194810, 0.399244},
    {61210, 0.0213557, 0.195458, 0.398585},
    {61211, 0.0213668, 0.196097, 0.397913},
    {61212, 0.0213030, 0.196728, 0.397230},
    {61213, 0.0212667, 0.197349, 0.396535},
    {61214, 0.0213450, 0.197959, 0.395829},
    {61215, 0.0216025, 0.198558, 0.395111},
    {61216, 0.0220766, 0.199145, 0.394382},
    {61217, 0.0227777, 0.199721, 0.393640},
    {61218, 0.0236872, 0.200284, 0.392886},
    {61219, 0.0247652, 0.200835, 0.392119},
    {61220, 0.0259447, 0.201374, 0.391342},
    {61221, 0.0271516, 0.201901, 0.390553},
    {61222, 0.0283109, 0.202414, 0.389752},
    {61223, 0.0293553, 0.202915, 0.388941},
    {61224, 0.0302303, 0.203402, 0.388120},
    {61225, 0.0309074, 0.203875, 0.387288},
    {61226, 0.0313892, 0.204335, 0.386446},
    {61227, 0.0317095, 0.204779, 0.385594},
    {61228, 0.0319410, 0.205209, 0.384732},
    {61229, 0.0321824, 0.205625, 0.383861},
    {61230, 0.0325458, 0.206024, 0.382980},
    {61231, 0.0331266, 0.206409, 0.382091},
    {61232, 0.0339635, 0.206777, 0.381192},
    {61233, 0.0350247, 0.207130, 0.380284},
    {61234, 0.0362043, 0.207466, 0.379368},
    {61235, 0.0373521, 0.207786, 0.378444},
    {61236, 0.0383277, 0.208090, 0.377511},
    {61237, 0.0390515, 0.208377, 0.376571},
    {61238, 0.0395386, 0.208646, 0.375624},
    {61239, 0.0398754, 0.208899, 0.374669},
    {61240, 0.0401837, 0.209133, 0.373707},
    {61241, 0.0405803, 0.209351, 0.372739},
    {61242, 0.0411450, 0.209550, 0.371764},
    {61243, 0.0419176, 0.209731, 0.370783},
    {61244, 0.0429050, 0.209894, 0.369796},
    {61245, 0.0440875, 0.210039, 0.368804},
    {61246, 0.0454240, 0.210165, 0.367806},
    {61247, 0.0468537, 0.210272, 0.366803},
    {61248, 0.0483024, 0.210361, 0.365795},
    {61249, 0.0496900, 0.210430, 0.364783},
    {61250, 0.0509413, 0.210480, 0.363767},
    {61251, 0.0519989, 0.210511, 0.362746},
    {61252, 0.0528299, 0.210523, 0.361722},
    {61253, 0.0534365, 0.210515, 0.360695},
    {61254, 0.0538612, 0.210487, 0.359665},
    {61255, 0.0541797, 0.210439, 0.358632},
    {61256, 0.0544894, 0.210371, 0.357597},
    {61257, 0.0548914, 0.210284, 0.356560},
    {61258, 0.0554700, 0.210176, 0.355521},
    {61259, 0.0562673, 0.210048, 0.354480},
    {61260, 0.0572647, 0.209900, 0.353439},
    {61261, 0.0583780, 0.209731, 0.352397},
    {61262, 0.0594771, 0.209542, 0.351354},
    {61263, 0.0604219, 0.209332, 0.350311},
    {61264, 0.0611086, 0.209102, 0.349268},
    {61265, 0.0615093, 0.208851, 0.348226},
    {61266, 0.0616755, 0.208580, 0.347184},
    {61267, 0.0617191, 0.208287, 0.346144},
    {61268, 0.0617732, 0.207975, 0.345105},
    {61269, 0.0619510, 0.207641, 0.344068},
    {61270, 0.0623211, 0.207287, 0.343033},
    {61271, 0.0629067, 0.206911, 0.342001},
    {61272, 0.0636950, 0.206515, 0.340971},
    {61273, 0.0646486, 0.206099, 0.339945},
    {61274, 0.0657097, 0.205661, 0.338922},
    {61275, 0.0668046, 0.205203, 0.337903},
    {61276, 0.0678550, 0.204724, 0.336888},
    {61277, 0.0687823, 0.204225, 0.335877},
    {61278, 0.0695182, 0.203704, 0.334872},
    {61279, 0.0700175, 0.203163, 0.333871},
    {61280, 0.0702675, 0.202602, 0.332876},
    {61281, 0.0702966, 0.202020, 0.331887},
    {61282, 0.0701776, 0.201418, 0.330904},
    {61283, 0.0700149, 0.200796, 0.329928},
    {61284, 0.0699228, 0.200155, 0.328959},
    {61285, 0.0699984, 0.199493, 0.327996},
    {61286, 0.0702939, 0.198811, 0.327042},
    {61287, 0.0707987, 0.198109, 0.326095},
    {61288, 0.0714425, 0.197388, 0.325156},
    {61289, 0.0721086, 0.196646, 0.324226},
    {61290, 0.0726304, 0.195885, 0.323305},
    {61291, 0.0729331, 0.195104, 0.322393},
    {61292, 0.0729585, 0.194304, 0.321490},
    {61293, 0.0727210, 0.193485, 0.320598},
    {61294, 0.0723023, 0.192646, 0.319715},
    {61295, 0.0718249, 0.191789, 0.318843},
    {61296, 0.0714141, 0.190913, 0.317981},
    {61297, 0.0711660, 0.190018, 0.317131},
    {61298, 0.0711299, 0.189104, 0.316292},
    {61299, 0.0713103, 0.188173, 0.315465},
    {61300, 0.0716767, 0.187223, 0.314650},
    {61301, 0.0721748, 0.186255, 0.313847},
    {61302, 0.0727357, 0.185270, 0.313057},
    {61303, 0.0732823, 0.184267, 0.312280},
    {61304, 0.0737362, 0.183247, 0.311516},
    {61305, 0.0740267, 0.182210, 0.310765},
    {61306, 0.0741005, 0.181156, 0.310029},
    {61307, 0.0739298, 0.180086, 0.309307},
    {61308, 0.0735245, 0.178999, 0.308599},
    {61309, 0.0729386, 0.177896, 0.307905},
    {61310, 0.0722666, 0.176777, 0.307227},
    {61311, 0.0716300, 0.175642, 0.306564},
    {61312, 0.0711462, 0.174492, 0.305917},
    {61313, 0.0708931, 0.173327, 0.305285},
    {61314, 0.0708821, 0.172147, 0.304670},
    {61315, 0.0710531, 0.170953, 0.304071},
    {61316, 0.0712938, 0.169744, 0.303488},
    {61317, 0.0714721, 0.168521, 0.302923},
    {61318, 0.0714739, 0.167285, 0.302374},
    {61319, 0.0712324, 0.166035, 0.301843},
    {61320, 0.0707421, 0.164772, 0.301330},
    {61321, 0.0700576, 0.163497, 0.300834},
    {61322, 0.0692781, 0.162208, 0.300356},
    {61323, 0.0685195, 0.160908, 0.299897},
    {61324, 0.0678860, 0.159596, 0.299456},
    {61325, 0.0674476, 0.158273, 0.299034},
    {61326, 0.0672311, 0.156938, 0.298631},
    {61327, 0.0672225, 0.155593, 0.298248},
    {61328, 0.0673768, 0.154237, 0.297883},
    {61329, 0.0676285, 0.152871, 0.297538},
    {61330, 0.0679029, 0.151495, 0.297213},
    {61331, 0.0681224, 0.150109, 0.296908},
    {61332, 0.0682160, 0.148715, 0.296623},
    {61333, 0.0681266, 0.147312, 0.296359},
    {61334, 0.0678172, 0.145901, 0.296115},
    {61335, 0.0672805, 0.144482, 0.295891},
    {61336, 0.0665454, 0.143055, 0.295689},
    {61337, 0.0656836, 0.141621, 0.295507},
    {61338, 0.0648042, 0.140180, 0.295346},
    {61339, 0.0640325, 0.138733, 0.295207},
    {61340, 0.0634737, 0.137280, 0.295090},
    {61341, 0.0631743, 0.135821, 0.294993},
    {61342, 0.0631014, 0.134357, 0.294919},
    {61343, 0.0631516, 0.132888, 0.294866},
    {61344, 0.0631864, 0.131414, 0.294835},
    {61345, 0.0630801, 0.129937, 0.294827},
    {61346, 0.0627561, 0.128456, 0.294840},
    {61347, 0.0622031, 0.126971, 0.294876},
    {61348, 0.0614679, 0.125484, 0.294933},
    {61349, 0.0606377, 0.123995, 0.295014},
    {61350, 0.0598146, 0.122504, 0.295117},
    {61351, 0.0590946, 0.121011, 0.295242},
    {61352, 0.0585497, 0.119517, 0.295390},
    {61353, 0.0582181, 0.118022, 0.295560},
    {61354, 0.0581007, 0.116527, 0.295754},
    {61355, 0.0581645, 0.115032, 0.295970},
    {61356, 0.0583517, 0.113538, 0.296208},
    {61357, 0.0585893, 0.112044, 0.296470},
    {61358, 0.0588003, 0.110552, 0.296754},
    {61359, 0.0589126, 0.109062, 0.297062},
    {61360, 0.0588670, 0.107575, 0.297392},
    {61361, 0.0586227, 0.106090, 0.297745},
    {61362, 0.0581617, 0.104608, 0.298121},
    {61363, 0.0574940, 0.103129, 0.298520},
    {61364, 0.0566625, 0.101655, 0.298941},
    {61365, 0.0557463, 0.100185, 0.299386},
    {61366, 0.0548514, 0.098720, 0.299853},
    {61367, 0.0540859, 0.097261, 0.300343},
    {61368, 0.0535234, 0.095807, 0.300856},
    {61369, 0.0531691, 0.094360, 0.301391},
    {61370, 0.0529502, 0.092919, 0.301949},
    {61371, 0.0527402, 0.091485, 0.302529},
    {61372, 0.0524081, 0.090058, 0.303132},
    {61373, 0.0518692, 0.088640, 0.303757},
    {61374, 0.0511116, 0.087230, 0.304405},
    {61375, 0.0501919, 0.085828, 0.305074},
    {61376, 0.0492078, 0.084436, 0.305766},
    {61377, 0.0482646, 0.083054, 0.306480},
    {61378, 0.0474512, 0.081681, 0.307215},
    {61379, 0.0468277, 0.080319, 0.307973},
    {61380, 0.0464240, 0.078968, 0.308752},
    {61381, 0.0462424, 0.077628, 0.309552},
    {61382, 0.0462611, 0.076300, 0.310374},
    {61383, 0.0464372, 0.074984, 0.311217},
    {61384, 0.0467111, 0.073681, 0.312081},
    {61385, 0.0470123, 0.072391, 0.312965},
    {61386, 0.0472680, 0.071114, 0.313871},
    {61387, 0.0474126, 0.069851, 0.314797},
    {61388, 0.0473981, 0.068602, 0.315743},
    {61389, 0.0472018, 0.067367, 0.316710},
    {61390, 0.0468310, 0.066148, 0.317697},
    {61391, 0.0463243, 0.064944, 0.318703},
    {61392, 0.0457495, 0.063756, 0.319729},
    {61393, 0.0451976, 0.062584, 0.320774},
    {61394, 0.0447684, 0.061428, 0.321838},
    {61395, 0.0445436, 0.060290, 0.322921},
    {61396, 0.0445541, 0.059168, 0.324023},
    {61397, 0.0447547, 0.058065, 0.325143},
    {61398, 0.0450225, 0.056979, 0.326282},
    {61399, 0.0451988, 0.055912, 0.327438},
    {61400, 0.0451697, 0.054863, 0.328612},
    {61401, 0.0449013, 0.053834, 0.329803},
    {61402, 0.0444463, 0.052824, 0.331012},
    {61403, 0.0439028, 0.051834, 0.332237},
    {61404, 0.0433500, 0.050863, 0.333479},
    {61405, 0.0428714, 0.049914, 0.334737},
};


static od_fit_vec3_t od_fit_vec_add(od_fit_vec3_t a, od_fit_vec3_t b)
{
    od_fit_vec3_t r;
    r.x = a.x + b.x;
    r.y = a.y + b.y;
    r.z = a.z + b.z;
    return r;
}

static double od_fit_anp(double angle)
{
    double r = fmod(angle, OD_FIT_D2PI);
    return (r < 0.0) ? r + OD_FIT_D2PI : r;
}

static void od_fit_mat_identity(double r[3][3])
{
    int i;
    int j;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            r[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

static void od_fit_mat_mul(const double a[3][3], const double b[3][3], double out[3][3])
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

static od_fit_vec3_t od_fit_mat_vec(const double m[3][3], od_fit_vec3_t v)
{
    od_fit_vec3_t r;
    r.x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
    r.y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
    r.z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;
    return r;
}

static void od_fit_rx(double phi, double r[3][3])
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

static void od_fit_ry(double theta, double r[3][3])
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

static void od_fit_rz(double psi, double r[3][3])
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

static void od_fit_rz_vector_matrix(double angle, double r[3][3])
{
    double s = sin(angle);
    double c = cos(angle);
    r[0][0] = c;  r[0][1] = -s; r[0][2] = 0.0;
    r[1][0] = s;  r[1][1] = c;  r[1][2] = 0.0;
    r[2][0] = 0.0; r[2][1] = 0.0; r[2][2] = 1.0;
}

static double od_fit_mjd_from_jd(double jd)
{
    return jd - 2400000.5;
}

static od_fit_eop_value_t od_fit_eop_from_record(const od_fit_eop_record_t *record)
{
    od_fit_eop_value_t value;
    value.dut1_seconds = record->dut1_seconds;
    value.x_pole_rad = record->x_pole_arcsec * OD_FIT_ARCSEC2RAD;
    value.y_pole_rad = record->y_pole_arcsec * OD_FIT_ARCSEC2RAD;
    return value;
}

static od_fit_eop_value_t od_fit_interpolate_eop(
    const od_fit_eop_record_t *a,
    const od_fit_eop_record_t *b,
    double mjd_utc)
{
    double span = (double)(b->mjd - a->mjd);
    double f = span > 0.0 ? (mjd_utc - (double)a->mjd) / span : 0.0;
    od_fit_eop_value_t value;
    value.dut1_seconds = a->dut1_seconds + f * (b->dut1_seconds - a->dut1_seconds);
    value.x_pole_rad = (a->x_pole_arcsec + f * (b->x_pole_arcsec - a->x_pole_arcsec)) * OD_FIT_ARCSEC2RAD;
    value.y_pole_rad = (a->y_pole_arcsec + f * (b->y_pole_arcsec - a->y_pole_arcsec)) * OD_FIT_ARCSEC2RAD;
    return value;
}

static od_fit_eop_value_t od_fit_lookup_eop(double mjd_utc)
{
    od_fit_eop_value_t empty;
    empty.dut1_seconds = 0.0;
    empty.x_pole_rad = 0.0;
    empty.y_pole_rad = 0.0;
#if !OD_FIT_USE_EOP_TABLE
    (void)mjd_utc;
    return empty;
#else
    size_t count = sizeof(od_eop_2026) / sizeof(od_eop_2026[0]);
    size_t i;
    if (count == 0U) {
        return empty;
    }
    if (mjd_utc <= (double)od_eop_2026[0].mjd) {
        return od_fit_eop_from_record(&od_eop_2026[0]);
    }
    if (mjd_utc >= (double)od_eop_2026[count - 1U].mjd) {
        return od_fit_eop_from_record(&od_eop_2026[count - 1U]);
    }
    for (i = 1U; i < count; ++i) {
        if (mjd_utc <= (double)od_eop_2026[i].mjd) {
            return od_fit_interpolate_eop(&od_eop_2026[i - 1U], &od_eop_2026[i], mjd_utc);
        }
    }
    return empty;
#endif
}

static void od_fit_polar_motion_matrix(double xp, double yp, double out[3][3])
{
    double ry_mat[3][3];
    double rx_mat[3][3];
    od_fit_mat_identity(ry_mat);
    od_fit_mat_identity(rx_mat);
    od_fit_ry(-xp, ry_mat);
    od_fit_rx(-yp, rx_mat);
    od_fit_mat_mul(ry_mat, rx_mat, out);
}

typedef struct od_fit_nut_term_t {
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
} od_fit_nut_term_t;

static void od_fit_nut00b(double jd_tt, double *dpsi, double *deps)
{
    static const od_fit_nut_term_t x[] = {
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
    const double u2r = OD_FIT_ARCSEC2RAD / 1.0e7;
    const double dp_plan = -0.135 * OD_FIT_MAS2RAD;
    const double de_plan = 0.388 * OD_FIT_MAS2RAD;
    double t = (jd_tt - OD_FIT_JD_J2000) / OD_FIT_JULIAN_CENTURY;
    double el = fmod(485868.249036 + 1717915923.2178 * t, 1296000.0) * OD_FIT_ARCSEC2RAD;
    double elp = fmod(1287104.79305 + 129596581.0481 * t, 1296000.0) * OD_FIT_ARCSEC2RAD;
    double f = fmod(335779.526232 + 1739527262.8478 * t, 1296000.0) * OD_FIT_ARCSEC2RAD;
    double d = fmod(1072260.70369 + 1602961601.2090 * t, 1296000.0) * OD_FIT_ARCSEC2RAD;
    double om = fmod(450160.398036 - 6962890.5431 * t, 1296000.0) * OD_FIT_ARCSEC2RAD;
    double dp = 0.0;
    double de = 0.0;
    int i;
    for (i = (int)(sizeof(x) / sizeof(x[0])) - 1; i >= 0; --i) {
        double arg = fmod((double)x[i].nl * el + (double)x[i].nlp * elp +
                          (double)x[i].nf * f + (double)x[i].nd * d +
                          (double)x[i].nom * om, OD_FIT_D2PI);
        double sarg = sin(arg);
        double carg = cos(arg);
        dp += (x[i].ps + x[i].pst * t) * sarg + x[i].pc * carg;
        de += (x[i].ec + x[i].ect * t) * carg + x[i].es * sarg;
    }
    *dpsi = dp * u2r + dp_plan;
    *deps = de * u2r + de_plan;
}

static void od_fit_bias_precession_matrix(double jd_tt, double rbp[3][3])
{
    const double eps0 = 84381.448 * OD_FIT_ARCSEC2RAD;
    const double dpsibi = -0.041775 * OD_FIT_ARCSEC2RAD;
    const double depsbi = -0.0068192 * OD_FIT_ARCSEC2RAD;
    const double dra0 = -0.0146 * OD_FIT_ARCSEC2RAD;
    double t = (jd_tt - OD_FIT_JD_J2000) / OD_FIT_JULIAN_CENTURY;
    double psia77 = (5038.7784 + (-1.07259 + (-0.001147) * t) * t) * t * OD_FIT_ARCSEC2RAD;
    double oma77 = eps0 + ((0.05127 + (-0.007726) * t) * t) * t * OD_FIT_ARCSEC2RAD;
    double chia = (10.5526 + (-2.38064 + (-0.001125) * t) * t) * t * OD_FIT_ARCSEC2RAD;
    double dpsipr = (-0.29965 * OD_FIT_ARCSEC2RAD) * t;
    double depspr = (-0.02524 * OD_FIT_ARCSEC2RAD) * t;
    double rb[3][3];
    double rp[3][3];
    double psia = psia77 + dpsipr;
    double oma = oma77 + depspr;

    od_fit_mat_identity(rb);
    od_fit_rz(dra0, rb);
    od_fit_ry(dpsibi * sin(eps0), rb);
    od_fit_rx(-depsbi, rb);

    od_fit_mat_identity(rp);
    od_fit_rx(eps0, rp);
    od_fit_rz(-psia, rp);
    od_fit_rx(-oma, rp);
    od_fit_rz(chia, rp);

    od_fit_mat_mul(rp, rb, rbp);
}

static double od_fit_obl80(double jd_tt)
{
    double t = (jd_tt - OD_FIT_JD_J2000) / OD_FIT_JULIAN_CENTURY;
    return OD_FIT_ARCSEC2RAD * (84381.448 + (-46.8150 + (-0.00059 + 0.001813 * t) * t) * t);
}

static void od_fit_pnm00b(double jd_tt, double rbpn[3][3], double *dpsi_out, double *epsa_out)
{
    double dpsi;
    double deps;
    double epsa;
    double rbp[3][3];
    double rn[3][3];
    double dpsipr;
    double depspr;
    od_fit_nut00b(jd_tt, &dpsi, &deps);
    dpsipr = (-0.29965 * OD_FIT_ARCSEC2RAD) * ((jd_tt - OD_FIT_JD_J2000) / OD_FIT_JULIAN_CENTURY);
    depspr = (-0.02524 * OD_FIT_ARCSEC2RAD) * ((jd_tt - OD_FIT_JD_J2000) / OD_FIT_JULIAN_CENTURY);
    (void)dpsipr;
    epsa = od_fit_obl80(jd_tt) + depspr;
    od_fit_bias_precession_matrix(jd_tt, rbp);
    od_fit_mat_identity(rn);
    od_fit_rx(epsa, rn);
    od_fit_rz(-dpsi, rn);
    od_fit_rx(-(epsa + deps), rn);
    od_fit_mat_mul(rn, rbp, rbpn);
    if (dpsi_out) {
        *dpsi_out = dpsi;
    }
    if (epsa_out) {
        *epsa_out = epsa;
    }
}

static double od_fit_era00(double jd_ut1)
{
    double d1 = floor(jd_ut1 - 0.5) + 0.5;
    double d2 = jd_ut1 - d1;
    double t = d1 + (d2 - OD_FIT_JD_J2000);
    double f = fmod(d1, 1.0) + fmod(d2, 1.0);
    return od_fit_anp(OD_FIT_D2PI * (f + 0.7790572732640 + 0.00273781191135448 * t));
}

static double od_fit_gmst00(double jd_ut1, double jd_tt)
{
    double t = (jd_tt - OD_FIT_JD_J2000) / OD_FIT_JULIAN_CENTURY;
    return od_fit_anp(od_fit_era00(jd_ut1) +
                  (0.014506 + (4612.15739966 +
                  (1.39667721 + (-0.00009344 + 0.00001882 * t) * t) * t) * t) *
                  OD_FIT_ARCSEC2RAD);
}

static void od_fit_itrs_to_j2000_matrix(const od_fit_time_t *time_utc, double out[3][3])
{
    double jd_utc = time_utc->jd_utc;
    double jd_tt = jd_utc + OD_FIT_TT_MINUS_UTC_SECONDS / OD_FIT_SECONDS_PER_DAY;
    od_fit_eop_value_t eop = od_fit_lookup_eop(od_fit_mjd_from_jd(jd_utc));
    double jd_ut1 = jd_utc + eop.dut1_seconds / OD_FIT_SECONDS_PER_DAY;
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
    od_fit_pnm00b(jd_tt, rbpn, &dpsi, &epsa);
    gst = od_fit_anp(od_fit_gmst00(jd_ut1, jd_tt) + dpsi * cos(epsa));
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            rbpn_t[i][j] = rbpn[j][i];
        }
    }
    od_fit_rz_vector_matrix(gst, rz);
    od_fit_polar_motion_matrix(-eop.x_pole_rad, -eop.y_pole_rad, pm);
    od_fit_mat_mul(rz, pm, rot_pm);
    od_fit_mat_mul(rbpn_t, rot_pm, out);
}

static void od_fit_ecef_to_j2000(
    const od_fit_time_t *time_utc,
    od_fit_vec3_t r_ecef,
    od_fit_vec3_t v_ecef,
    od_fit_vec3_t *r_j2000,
    od_fit_vec3_t *v_j2000)
{
    double m[3][3];
    od_fit_vec3_t spin;
    od_fit_vec3_t v_total;
    od_fit_itrs_to_j2000_matrix(time_utc, m);
    spin.x = -OD_FIT_OMEGA_EARTH * r_ecef.y;
    spin.y = OD_FIT_OMEGA_EARTH * r_ecef.x;
    spin.z = 0.0;
    v_total = od_fit_vec_add(v_ecef, spin);
    *r_j2000 = od_fit_mat_vec(m, r_ecef);
    *v_j2000 = od_fit_mat_vec(m, v_total);
}

static od_fit_vec3_t od_fit_observation_to_j2000_position(const od_fit_observation_t *observation)
{
    od_fit_vec3_t zero;
    od_fit_vec3_t r_j2000;
    od_fit_vec3_t ignored_v_j2000;
    zero.x = zero.y = zero.z = 0.0;
    od_fit_ecef_to_j2000(
        &observation->time_utc,
        observation->r_ecef_m,
        zero,
        &r_j2000,
        &ignored_v_j2000);
    return r_j2000;
}


static od_fit_options_t od_fit_default_options(void)
{
    od_fit_options_t opt;
    opt.degree = 10;
    return opt;
}

static void od_fit_cheb_basis(double tau, int degree, double *basis)
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

static int od_fit_solve_linear(int n, double a[OD_FIT_MAX_DEGREE + 1][OD_FIT_MAX_DEGREE + 2], double *x)
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
            return OD_FIT_ERR_FIT;
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
    return OD_FIT_OK;
}

static int od_fit_build_fit(
    const od_fit_observation_t *obs,
    size_t first,
    size_t last,
    int requested_degree,
    od_fit_fit_t *fit)
{
    size_t nobs;
    int degree;
    int ncoef;
    double ata[OD_FIT_MAX_DEGREE + 1][OD_FIT_MAX_DEGREE + 1];
    double rhs[3][OD_FIT_MAX_DEGREE + 1];
    size_t i;
    int j;
    int k;
    int coord;
    double basis[OD_FIT_MAX_DEGREE + 1];
    if (!obs || !fit || last < first) {
        return OD_FIT_ERR_INVALID_ARGUMENT;
    }
    nobs = last - first + 1;
    if (nobs < 2) {
        return OD_FIT_ERR_FIT;
    }
    degree = requested_degree;
    if (degree < 1) {
        degree = 1;
    }
    if (degree > OD_FIT_MAX_DEGREE) {
        degree = OD_FIT_MAX_DEGREE;
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
        return OD_FIT_ERR_FIT;
    }
    for (i = first; i <= last; ++i) {
        double tau = (obs[i].time_utc.unix_seconds - fit->epoch_seconds) / fit->scale_seconds;
        od_fit_vec3_t r_j2000 = od_fit_observation_to_j2000_position(&obs[i]);
        od_fit_cheb_basis(tau, degree, basis);
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
        double aug[OD_FIT_MAX_DEGREE + 1][OD_FIT_MAX_DEGREE + 2];
        double x[OD_FIT_MAX_DEGREE + 1];
        int status;
        for (j = 0; j < ncoef; ++j) {
            for (k = 0; k < ncoef; ++k) {
                aug[j][k] = ata[j][k];
            }
            aug[j][ncoef] = rhs[coord][j];
            x[j] = 0.0;
        }
        status = od_fit_solve_linear(ncoef, aug, x);
        if (status != OD_FIT_OK) {
            return status;
        }
        for (j = 0; j < ncoef; ++j) {
            fit->coeff[coord][j] = x[j];
        }
    }
    return OD_FIT_OK;
}

static double od_fit_cheb_eval(const double *coeff, int degree, double tau)
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

static double od_fit_cheb_derivative_eval(const double *coeff, int degree, double tau)
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

static void od_fit_fit_eval(const od_fit_fit_t *fit, const od_fit_time_t *query, od_fit_state_t *out)
{
    double tau = (query->unix_seconds - fit->epoch_seconds) / fit->scale_seconds;
    out->time_utc = *query;
    out->r_j2000_m.x = od_fit_cheb_eval(fit->coeff[0], fit->degree, tau);
    out->r_j2000_m.y = od_fit_cheb_eval(fit->coeff[1], fit->degree, tau);
    out->r_j2000_m.z = od_fit_cheb_eval(fit->coeff[2], fit->degree, tau);
    out->v_j2000_mps.x = od_fit_cheb_derivative_eval(fit->coeff[0], fit->degree, tau) / fit->scale_seconds;
    out->v_j2000_mps.y = od_fit_cheb_derivative_eval(fit->coeff[1], fit->degree, tau) / fit->scale_seconds;
    out->v_j2000_mps.z = od_fit_cheb_derivative_eval(fit->coeff[2], fit->degree, tau) / fit->scale_seconds;
}

static int od_fit_interpolate_state_cached(
    const od_fit_observation_t *obs,
    size_t count,
    const od_fit_time_t *query,
    const od_fit_options_t *opt,
    od_fit_fit_cache_t *cache,
    od_fit_state_t *out)
{
    size_t first;
    size_t last;
    int degree = opt->degree;
    int status;
    first = 0;
    last = count - 1;
    if (!cache || !cache->valid || cache->first_index != first || cache->last_index != last || cache->degree != degree) {
        od_fit_fit_t fit;
        status = od_fit_build_fit(obs, first, last, degree, &fit);
        if (status != OD_FIT_OK) {
            return status;
        }
        if (cache) {
            cache->valid = 1;
            cache->degree = degree;
            cache->first_index = first;
            cache->last_index = last;
            cache->fit = fit;
        } else {
            od_fit_fit_eval(&fit, query, out);
            return OD_FIT_OK;
        }
    }
    od_fit_fit_eval(&cache->fit, query, out);
    return OD_FIT_OK;
}

static void od_fit_context_invalidate_cache(od_fit_context_t *context)
{
    memset(&context->main_cache, 0, sizeof(context->main_cache));
}

static void od_fit_reverse_observations(od_fit_observation_t *obs, size_t first, size_t last)
{
    while (first < last) {
        od_fit_observation_t tmp = obs[first];
        obs[first] = obs[last];
        obs[last] = tmp;
        ++first;
        --last;
    }
}

static void od_fit_context_linearize(od_fit_context_t *context)
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

    od_fit_reverse_observations(context->observations, 0, split - 1);
    od_fit_reverse_observations(context->observations, split, context->count - 1);
    od_fit_reverse_observations(context->observations, 0, context->count - 1);
    context->start = 0;
    od_fit_context_invalidate_cache(context);
}

static int od_fit_context_create(
    od_fit_context_t **out_context,
    od_fit_observation_t *observation_buffer,
    size_t capacity,
    const od_fit_options_t *options)
{
    od_fit_context_t *context;
    if (!out_context || !observation_buffer) {
        return OD_FIT_ERR_INVALID_ARGUMENT;
    }
    *out_context = NULL;
    if (capacity == 0) {
        capacity = OD_FIT_DEFAULT_OBSERVATION_CAPACITY;
    }
    if (capacity < 2) {
        return OD_FIT_ERR_INVALID_ARGUMENT;
    }

    context = (od_fit_context_t *)malloc(sizeof(*context));
    if (!context) {
        return OD_FIT_ERR_NO_MEMORY;
    }
    memset(context, 0, sizeof(*context));
    context->observations = observation_buffer;
    context->capacity = capacity;
    context->options = options ? *options : od_fit_default_options();
    od_fit_context_invalidate_cache(context);
    *out_context = context;
    return OD_FIT_OK;
}

static void od_fit_context_reset(od_fit_context_t *context)
{
    if (!context) {
        return;
    }
    context->count = 0;
    context->start = 0;
    od_fit_context_invalidate_cache(context);
}

static void od_fit_context_destroy(od_fit_context_t *context)
{
    if (!context) {
        return;
    }
    memset(context, 0, sizeof(*context));
    free(context);
}

static size_t od_fit_context_count(const od_fit_context_t *context)
{
    return context ? context->count : 0;
}


static int od_fit_context_push(
    od_fit_context_t *context,
    const od_fit_observation_t *observation)
{
    size_t index;
    if (!context || !context->observations || context->capacity < 2 || !observation ||
        !isfinite(observation->r_ecef_m.x) ||
        !isfinite(observation->r_ecef_m.y) ||
        !isfinite(observation->r_ecef_m.z)) {
        return OD_FIT_ERR_INVALID_ARGUMENT;
    }
    if (context->count > 0) {
        size_t latest_index = (context->start + context->count - 1) % context->capacity;
        double latest_time = context->observations[latest_index].time_utc.unix_seconds;
        if (observation->time_utc.unix_seconds <= latest_time) {
            return OD_FIT_ERR_RANGE;
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

    od_fit_context_invalidate_cache(context);
    return OD_FIT_OK;
}

static int od_fit_context_query_state(
    od_fit_context_t *context,
    const od_fit_time_t *query_time_utc,
    od_fit_state_t *out_state)
{
    if (!context || !context->observations || !query_time_utc || !out_state) {
        return OD_FIT_ERR_INVALID_ARGUMENT;
    }
    if (context->count < 2) {
        return OD_FIT_ERR_INVALID_ARGUMENT;
    }

    od_fit_context_linearize(context);
    if (query_time_utc->unix_seconds < context->observations[0].time_utc.unix_seconds ||
        query_time_utc->unix_seconds > context->observations[context->count - 1].time_utc.unix_seconds) {
        return OD_FIT_ERR_RANGE;
    }
    return od_fit_interpolate_state_cached(context->observations, context->count,
                                      query_time_utc, &context->options,
                                      &context->main_cache, out_state);
}


#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wdiscarded-qualifiers"
#endif
/*

The RTKLIB software package is distributed under the following BSD 2-clause
license. Users are permitted to develop, produce or sell their own non-
commercial or commercial products utilizing, linking or including RTKLIB as long
as they comply with the license.

--------------------------------------------------------------------------------

         Copyright (c) 2007-2020, T. Takasu, All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer. Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the
distribution.

The software package includes some companion executive binaries or shared
libraries necessary to execute APs on Windows. These licenses succeed to the
original ones of these software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

--------------------------------------------------------------------------------

Notes:
Previous versions of RTKLIB until ver. 2.4.1 had been distributed under GPLv3
license.


*/
/* Embedded RTKLIB 2.4.3 b34. Original copyright notices follow. */

#define NO_STRICT
#define ENAGLO
#define ENAQZS
#define ENAGAL
#define ENACMP
#define ENAIRN
#define NFREQ 5
#define NEXOBS 3

/*------------------------------------------------------------------------------
* rtklib.h : RTKLIB constants, types and function prototypes
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* options : -DENAGLO   enable GLONASS
*           -DENAGAL   enable Galileo
*           -DENAQZS   enable QZSS
*           -DENACMP   enable BeiDou
*           -DENAIRN   enable IRNSS
*           -DNFREQ=n  set number of obs codes/frequencies
*           -DNEXOBS=n set number of extended obs codes
*           -DMAXOBS=n set max number of obs data in an epoch
*           -DWIN32    use WIN32 API
*           -DWIN_DLL  generate library as Windows DLL
*
* version : $Revision:$ $Date:$
* history : 2007/01/13 1.0  rtklib ver.1.0.0
*           2007/03/20 1.1  rtklib ver.1.1.0
*           2008/07/15 1.2  rtklib ver.2.1.0
*           2008/10/19 1.3  rtklib ver.2.1.1
*           2009/01/31 1.4  rtklib ver.2.2.0
*           2009/04/30 1.5  rtklib ver.2.2.1
*           2009/07/30 1.6  rtklib ver.2.2.2
*           2009/12/25 1.7  rtklib ver.2.3.0
*           2010/07/29 1.8  rtklib ver.2.4.0
*           2011/05/27 1.9  rtklib ver.2.4.1
*           2013/03/28 1.10 rtklib ver.2.4.2
*           2020/11/30 1.11 rtklib ver.2.4.3 b34
*-----------------------------------------------------------------------------*/
#ifndef RTKLIB_H
#define RTKLIB_H
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <stdint.h>
#ifdef WIN32
#include <winsock2.h>
#include <windows.h>
#else
#include <pthread.h>
#include <sys/select.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN_DLL
#define EXPORT __declspec(dllexport) /* for Windows DLL */
#else
#define EXPORT
#endif

/* constants -----------------------------------------------------------------*/

#define VER_RTKLIB  "2.4.3"             /* library version */

#define PATCH_LEVEL "b34"               /* patch level */

#define COPYRIGHT_RTKLIB \
            "Copyright (C) 2007-2020 T.Takasu\nAll rights reserved."

#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */

#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */

#define HION        350000.0            /* ionosphere height (m) */

#define MAXFREQ     7                   /* max NFREQ */

#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2         frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a/B2a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/L6  frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9          /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ1a_GLO  1.600995E9          /* GLONASS G1a frequency (Hz) */
#define FREQ2a_GLO  1.248060E9          /* GLONASS G2a frequency (Hz) */
#define FREQ1_CMP   1.561098E9          /* BDS B1I     frequency (Hz) */
#define FREQ2_CMP   1.20714E9           /* BDS B2I/B2b frequency (Hz) */
#define FREQ3_CMP   1.26852E9           /* BDS B3      frequency (Hz) */

#define EFACT_GPS   1.0                 /* error factor: GPS */
#define EFACT_GLO   1.5                 /* error factor: GLONASS */
#define EFACT_GAL   1.0                 /* error factor: Galileo */
#define EFACT_QZS   1.0                 /* error factor: QZSS */
#define EFACT_CMP   1.0                 /* error factor: BeiDou */
#define EFACT_IRN   1.5                 /* error factor: IRNSS */
#define EFACT_SBS   3.0                 /* error factor: SBAS */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_CMP     0x20                /* navigation system: BeiDou */
#define SYS_IRN     0x40                /* navigation system: IRNS */
#define SYS_LEO     0x80                /* navigation system: LEO */
#define SYS_ALL     0xFF                /* navigation system: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
#define TSYS_IRN    6                   /* time system: IRNSS time */

#ifndef NFREQ
#define NFREQ       3                   /* number of carrier frequencies */
#endif
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */

#ifndef NEXOBS
#define NEXOBS      0                   /* number of extended obs codes */
#endif

#define SNR_UNIT    0.001               /* SNR unit (dBHz) */

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#ifdef ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif
#ifdef ENAGAL
#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   36                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1
#else
#define MINPRNGAL   0
#define MAXPRNGAL   0
#define NSATGAL     0
#define NSYSGAL     0
#endif
#ifdef ENAQZS
#define MINPRNQZS   193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS   202                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                 /* min satellite PRN number of QZSS L1S */
#define MAXPRNQZS_S 191                 /* max satellite PRN number of QZSS L1S */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif
#ifdef ENACMP
#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   63                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1
#else
#define MINPRNCMP   0
#define MAXPRNCMP   0
#define NSATCMP     0
#define NSYSCMP     0
#endif
#ifdef ENAIRN
#define MINPRNIRN   1                   /* min satellite sat number of IRNSS */
#define MAXPRNIRN   14                  /* max satellite sat number of IRNSS */
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1) /* number of IRNSS satellites */
#define NSYSIRN     1
#else
#define MINPRNIRN   0
#define MAXPRNIRN   0
#define NSATIRN     0
#define NSYSIRN     0
#endif
#ifdef ENALEO
#define MINPRNLEO   1                   /* min satellite sat number of LEO */
#define MAXPRNLEO   10                  /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif
#define NSYS        (NSYSGPS+NSYSGLO+NSYSGAL+NSYSQZS+NSYSCMP+NSYSIRN+NSYSLEO) /* number of systems */

#define MINPRNSBS   120                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   158                 /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */

#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
                                        /* max satellite number (1 to MAXSAT) */
#define MAXSTA      255

#ifndef MAXOBS
#define MAXOBS      96                  /* max number of obs in an epoch */
#endif
#define MAXRCV      64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */
#ifdef OBS_100HZ
#define DTTOL       0.005               /* tolerance of time difference (s) */
#else
#define DTTOL       0.025               /* tolerance of time difference (s) */
#endif
#define MAXDTOE     7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_IRN 7200.0              /* max time difference to IRNSS Toe (s) */
#define MAXDTOE_SBS 360.0               /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S   86400.0             /* max time difference to ephem toe (s) for other */
#define MAXGDOP     300.0               /* max GDOP */

#define INT_SWAP_TRAC 86400.0           /* swap interval of trace file (s) */
#define INT_SWAP_STAT 86400.0           /* swap interval of solution status file (s) */

#define MAXEXFILE   1024                /* max number of expanded files */
#define MAXSBSAGEF  30.0                /* max age of SBAS fast correction (s) */
#define MAXSBSAGEL  1800.0              /* max age of SBAS long term corr (s) */
#define MAXSBSURA   8                   /* max URA of SBAS satellite */
#define MAXBAND     10                  /* max SBAS band of IGP */
#define MAXNIGP     201                 /* max number of IGP in SBAS band */
#define MAXNGEO     4                   /* max number of GEO satellites */
#define MAXCOMMENT  100                 /* max number of RINEX comments */
#define MAXSTRPATH  1024                /* max length of stream path */
#define MAXSTRMSG   1024                /* max length of stream message */
#define MAXSTRRTK   8                   /* max number of stream in RTK server */
#define MAXSBSMSG   32                  /* max number of SBAS msg in RTK server */
#define MAXSOLMSG   8191                /* max length of solution message */
#define MAXRAWLEN   16384               /* max length of receiver raw message */
#define MAXERRMSG   4096                /* max length of error/warning message */
#define MAXANT      64                  /* max length of station name/antenna type */
#define MAXSOLBUF   256                 /* max number of solution buffer */
#define MAXOBSBUF   128                 /* max number of observation data buffer */
#define MAXNRPOS    16                  /* max number of reference positions */
#define MAXLEAPS    64                  /* max number of leap seconds table */
#define MAXGISLAYER 32                  /* max number of GIS data layers */
#define MAXRCVCMD   4096                /* max length of receiver commands */

#define RNX2VER     2.10                /* RINEX ver.2 default output version */
#define RNX3VER     3.00                /* RINEX ver.3 default output version */

#define OBSTYPE_PR  0x01                /* observation type: pseudorange */
#define OBSTYPE_CP  0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP 0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR 0x08                /* observation type: SNR */
#define OBSTYPE_ALL 0xFF                /* observation type: all */

#define FREQTYPE_L1 0x01                /* frequency type: L1/E1/B1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/E5b/B2 */
#define FREQTYPE_L3 0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_L4 0x08                /* frequency type: L6/E6/B3 */
#define FREQTYPE_L5 0x10                /* frequency type: E5ab */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless,B1codeless (GPS,BDS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A,B1A    (GAL,BDS) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1S (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2bI  (GAL,BDS) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2bQ  (GAL,BDS) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
#define CODE_L6A    30                  /* obs code: E6A,B3A    (GAL,BDS) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
#define CODE_L6S    35                  /* obs code: L6S        (QZS) */
#define CODE_L6L    36                  /* obs code: L6L        (QZS) */
#define CODE_L8I    37                  /* obs code: E5abI      (GAL) */
#define CODE_L8Q    38                  /* obs code: E5abQ      (GAL) */
#define CODE_L8X    39                  /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
#define CODE_L2I    40                  /* obs code: B1_2I      (BDS) */
#define CODE_L2Q    41                  /* obs code: B1_2Q      (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) (obsolute) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) (obsolute) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define CODE_L1D    56                  /* obs code: B1D        (BDS) */
#define CODE_L5D    57                  /* obs code: L5D(L5S),B2aD (QZS,BDS) */
#define CODE_L5P    58                  /* obs code: L5P(L5S),B2aP (QZS,BDS) */
#define CODE_L5Z    59                  /* obs code: L5D+P(L5S) (QZS) */
#define CODE_L6E    60                  /* obs code: L6E        (QZS) */
#define CODE_L7D    61                  /* obs code: B2bD       (BDS) */
#define CODE_L7P    62                  /* obs code: B2bP       (BDS) */
#define CODE_L7Z    63                  /* obs code: B2bD+P     (BDS) */
#define CODE_L8D    64                  /* obs code: B2abD      (BDS) */
#define CODE_L8P    65                  /* obs code: B2abP      (BDS) */
#define CODE_L4A    66                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4B    67                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4X    68                  /* obs code: G1al1OCd+p (GLO) */
#define MAXCODE     68                  /* max number of obs code */

#define PMODE_SINGLE 0                  /* positioning mode: single */
#define PMODE_DGPS   1                  /* positioning mode: DGPS/DGNSS */
#define PMODE_KINEMA 2                  /* positioning mode: kinematic */
#define PMODE_STATIC 3                  /* positioning mode: static */
#define PMODE_MOVEB  4                  /* positioning mode: moving-base */
#define PMODE_FIXED  5                  /* positioning mode: fixed */
#define PMODE_PPP_KINEMA 6              /* positioning mode: PPP-kinemaric */
#define PMODE_PPP_STATIC 7              /* positioning mode: PPP-static */
#define PMODE_PPP_FIXED 8               /* positioning mode: PPP-fixed */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */
#define SOLF_ENU    2                   /* solution format: e/n/u-baseline */
#define SOLF_NMEA   3                   /* solution format: NMEA-183 */
#define SOLF_STAT   4                   /* solution format: solution status */
#define SOLF_GSIF   5                   /* solution format: GSI F1/F2 */

#define SOLQ_NONE   0                   /* solution status: no solution */
#define SOLQ_FIX    1                   /* solution status: fix */
#define SOLQ_FLOAT  2                   /* solution status: float */
#define SOLQ_SBAS   3                   /* solution status: SBAS */
#define SOLQ_DGPS   4                   /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 5                   /* solution status: single */
#define SOLQ_PPP    6                   /* solution status: PPP */
#define SOLQ_DR     7                   /* solution status: dead reconing */
#define MAXSOLQ     7                   /* max number of solution status */

#define TIMES_GPST  0                   /* time system: gps time */
#define TIMES_UTC   1                   /* time system: utc */
#define TIMES_JST   2                   /* time system: jst */

#define IONOOPT_OFF 0                   /* ionosphere option: correction off */
#define IONOOPT_BRDC 1                  /* ionosphere option: broadcast model */
#define IONOOPT_SBAS 2                  /* ionosphere option: SBAS model */
#define IONOOPT_IFLC 3                  /* ionosphere option: L1/L2 iono-free LC */
#define IONOOPT_EST 4                   /* ionosphere option: estimation */
#define IONOOPT_TEC 5                   /* ionosphere option: IONEX TEC model */
#define IONOOPT_QZS 6                   /* ionosphere option: QZSS broadcast model */
#define IONOOPT_STEC 8                  /* ionosphere option: SLANT TEC model */

#define TROPOPT_OFF 0                   /* troposphere option: correction off */
#define TROPOPT_SAAS 1                  /* troposphere option: Saastamoinen model */
#define TROPOPT_SBAS 2                  /* troposphere option: SBAS model */
#define TROPOPT_EST 3                   /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG 4                  /* troposphere option: ZTD+grad estimation */
#define TROPOPT_ZTD 5                   /* troposphere option: ZTD correction */

#define EPHOPT_BRDC 0                   /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC 1                   /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS 2                   /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC 3                 /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM 4                 /* ephemeris option: broadcast + SSR_COM */

#define ARMODE_OFF  0                   /* AR mode: off */
#define ARMODE_CONT 1                   /* AR mode: continuous */
#define ARMODE_INST 2                   /* AR mode: instantaneous */
#define ARMODE_FIXHOLD 3                /* AR mode: fix and hold */
#define ARMODE_WLNL 4                   /* AR mode: wide lane/narrow lane */
#define ARMODE_TCAR 5                   /* AR mode: triple carrier ar */

#define SBSOPT_LCORR 1                  /* SBAS option: long term correction */
#define SBSOPT_FCORR 2                  /* SBAS option: fast correction */
#define SBSOPT_ICORR 4                  /* SBAS option: ionosphere correction */
#define SBSOPT_RANGE 8                  /* SBAS option: ranging */

#define POSOPT_POS   0                  /* pos option: LLH/XYZ */
#define POSOPT_SINGLE 1                 /* pos option: average of single pos */
#define POSOPT_FILE  2                  /* pos option: read from pos file */
#define POSOPT_RINEX 3                  /* pos option: rinex header pos */
#define POSOPT_RTCM  4                  /* pos option: rtcm/raw station pos */

#define STR_NONE     0                  /* stream type: none */
#define STR_SERIAL   1                  /* stream type: serial */
#define STR_FILE     2                  /* stream type: file */
#define STR_TCPSVR   3                  /* stream type: TCP server */
#define STR_TCPCLI   4                  /* stream type: TCP client */
#define STR_NTRIPSVR 5                  /* stream type: NTRIP server */
#define STR_NTRIPCLI 6                  /* stream type: NTRIP client */
#define STR_FTP      7                  /* stream type: ftp */
#define STR_HTTP     8                  /* stream type: http */
#define STR_NTRIPCAS 9                  /* stream type: NTRIP caster */
#define STR_UDPSVR   10                 /* stream type: UDP server */
#define STR_UDPCLI   11                 /* stream type: UDP server */
#define STR_MEMBUF   12                 /* stream type: memory buffer */

#define STRFMT_RTCM2 0                  /* stream format: RTCM 2 */
#define STRFMT_RTCM3 1                  /* stream format: RTCM 3 */
#define STRFMT_OEM4  2                  /* stream format: NovAtel OEMV/4 */
#define STRFMT_OEM3  3                  /* stream format: NovAtel OEM3 */
#define STRFMT_UBX   4                  /* stream format: u-blox LEA-*T */
#define STRFMT_SS2   5                  /* stream format: NovAtel Superstar II */
#define STRFMT_CRES  6                  /* stream format: Hemisphere */
#define STRFMT_STQ   7                  /* stream format: SkyTraq S1315F */
#define STRFMT_JAVAD 8                  /* stream format: JAVAD GRIL/GREIS */
#define STRFMT_NVS   9                  /* stream format: NVS NVC08C */
#define STRFMT_BINEX 10                 /* stream format: BINEX */
#define STRFMT_RT17  11                 /* stream format: Trimble RT17 */
#define STRFMT_SEPT  12                 /* stream format: Septentrio */
#define STRFMT_RINEX 13                 /* stream format: RINEX */
#define STRFMT_SP3   14                 /* stream format: SP3 */
#define STRFMT_RNXCLK 15                /* stream format: RINEX CLK */
#define STRFMT_SBAS  16                 /* stream format: SBAS messages */
#define STRFMT_NMEA  17                 /* stream format: NMEA 0183 */
#define MAXRCVFMT    12                 /* max number of receiver format */

#define STR_MODE_R  0x1                 /* stream mode: read */
#define STR_MODE_W  0x2                 /* stream mode: write */
#define STR_MODE_RW 0x3                 /* stream mode: read/write */

#define GEOID_EMBEDDED    0             /* geoid model: embedded geoid */
#define GEOID_EGM96_M150  1             /* geoid model: EGM96 15x15" */
#define GEOID_EGM2008_M25 2             /* geoid model: EGM2008 2.5x2.5" */
#define GEOID_EGM2008_M10 3             /* geoid model: EGM2008 1.0x1.0" */
#define GEOID_GSI2000_M15 4             /* geoid model: GSI geoid 2000 1.0x1.5" */
#define GEOID_RAF09       5             /* geoid model: IGN RAF09 for France 1.5"x2" */

#define COMMENTH    "%"                 /* comment line indicator for solution */
#define MSG_DISCONN "$_DISCONNECT\r\n"  /* disconnect message */

#define DLOPT_FORCE   0x01              /* download option: force download existing */
#define DLOPT_KEEPCMP 0x02              /* download option: keep compressed file */
#define DLOPT_HOLDERR 0x04              /* download option: hold on error file */
#define DLOPT_HOLDLST 0x08              /* download option: hold on listing file */

#define LLI_SLIP    0x01                /* LLI: cycle-slip */
#define LLI_HALFC   0x02                /* LLI: half-cycle not resovled */
#define LLI_BOCTRK  0x04                /* LLI: boc tracking of mboc signal */
#define LLI_HALFA   0x40                /* LLI: half-cycle added */
#define LLI_HALFS   0x80                /* LLI: half-cycle subtracted */

#define P2_5        0.03125             /* 2^-5 */
#define P2_6        0.015625            /* 2^-6 */
#define P2_11       4.882812500000000E-04 /* 2^-11 */
#define P2_15       3.051757812500000E-05 /* 2^-15 */
#define P2_17       7.629394531250000E-06 /* 2^-17 */
#define P2_19       1.907348632812500E-06 /* 2^-19 */
#define P2_20       9.536743164062500E-07 /* 2^-20 */
#define P2_21       4.768371582031250E-07 /* 2^-21 */
#define P2_23       1.192092895507810E-07 /* 2^-23 */
#define P2_24       5.960464477539063E-08 /* 2^-24 */
#define P2_27       7.450580596923828E-09 /* 2^-27 */
#define P2_29       1.862645149230957E-09 /* 2^-29 */
#define P2_30       9.313225746154785E-10 /* 2^-30 */
#define P2_31       4.656612873077393E-10 /* 2^-31 */
#define P2_32       2.328306436538696E-10 /* 2^-32 */
#define P2_33       1.164153218269348E-10 /* 2^-33 */
#define P2_35       2.910383045673370E-11 /* 2^-35 */
#define P2_38       3.637978807091710E-12 /* 2^-38 */
#define P2_39       1.818989403545856E-12 /* 2^-39 */
#define P2_40       9.094947017729280E-13 /* 2^-40 */
#define P2_43       1.136868377216160E-13 /* 2^-43 */
#define P2_48       3.552713678800501E-15 /* 2^-48 */
#define P2_50       8.881784197001252E-16 /* 2^-50 */
#define P2_55       2.775557561562891E-17 /* 2^-55 */

#ifdef WIN32
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#define FILEPATHSEP '\\'
#else
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define lock(f)     pthread_mutex_lock(f)
#define unlock(f)   pthread_mutex_unlock(f)
#define FILEPATHSEP '/'
#endif

/* type definitions ----------------------------------------------------------*/

typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
} gtime_t;

typedef struct {        /* observation data record */
    gtime_t time;       /* receiver sampling time (GPST) */
    uint8_t sat,rcv;    /* satellite/receiver number */
    uint16_t SNR[NFREQ+NEXOBS]; /* signal strength (0.001 dBHz) */
    uint8_t  LLI[NFREQ+NEXOBS]; /* loss of lock indicator */
    uint8_t code[NFREQ+NEXOBS]; /* code indicator (CODE_???) */
    double L[NFREQ+NEXOBS]; /* observation data carrier-phase (cycle) */
    double P[NFREQ+NEXOBS]; /* observation data pseudorange (m) */
    float  D[NFREQ+NEXOBS]; /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {        /* observation data */
    int n,nmax;         /* number of obervation data/allocated */
    obsd_t *data;       /* observation data records */
} obs_t;

typedef struct {        /* earth rotation parameter data type */
    double mjd;         /* mjd (days) */
    double xp,yp;       /* pole offset (rad) */
    double xpr,ypr;     /* pole offset rate (rad/day) */
    double ut1_utc;     /* ut1-utc (s) */
    double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
    int n,nmax;         /* number and max number of data */
    erpd_t *data;       /* earth rotation parameter data */
} erp_t;

typedef struct {        /* antenna parameter type */
    int sat;            /* satellite number (0:receiver) */
    char type[MAXANT];  /* antenna type */
    char code[MAXANT];  /* serial number or satellite code */
    gtime_t ts,te;      /* valid time start and end */
    double off[NFREQ][ 3]; /* phase center offset e/n/u or x/y/z (m) */
    double var[NFREQ][19]; /* phase center variation (m) */
                        /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} pcv_t;

typedef struct {        /* antenna parameters type */
    int n,nmax;         /* number of data/allocated */
    pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

typedef struct {        /* almanac type */
    int sat;            /* satellite number */
    int svh;            /* sv health (0:ok) */
    int svconf;         /* as and sv config */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    gtime_t toa;        /* Toa */
                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,OMGd;
    double toas;        /* Toa (s) in week */
    double f0,f1;       /* SV clock parameters (af0,af1) */
} alm_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode,iodc;      /* IODE,IODC */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    int code;           /* GPS/QZS: code on L2 */
                        /* GAL: data source defined as rinex 3.03 */
                        /* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
    int flag;           /* GPS/QZS: L2 P data flag */
                        /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    double crc,crs,cuc,cus,cic,cis;
    double toes;        /* Toe (s) in week */
    double fit;         /* fit interval (h) */
    double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[6];      /* group delay parameters */
                        /* GPS/QZS:tgd[0]=TGD */
                        /* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
                        /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
                        /*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
    double Adot,ndot;   /* Adot,ndot for CNAV */
} eph_t;

typedef struct {        /* GLONASS broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode;           /* IODE (0-6 bit of tb field) */
    int frq;            /* satellite frequency number */
    int svh,sva,age;    /* satellite health, accuracy, age of operation */
    gtime_t toe;        /* epoch of epherides (gpst) */
    gtime_t tof;        /* message frame time (gpst) */
    double pos[3];      /* satellite position (ecef) (m) */
    double vel[3];      /* satellite velocity (ecef) (m/s) */
    double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;   /* SV clock bias (s)/relative freq bias */
    double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

typedef struct {        /* precise ephemeris type */
    gtime_t time;       /* time (GPST) */
    int index;          /* ephemeris index for multiple files */
    double pos[MAXSAT][4]; /* satellite position/clock (ecef) (m|s) */
    float  std[MAXSAT][4]; /* satellite position/clock std (m|s) */
    double vel[MAXSAT][4]; /* satellite velocity/clk-rate (m/s|s/s) */
    float  vst[MAXSAT][4]; /* satellite velocity/clk-rate std (m/s|s/s) */
    float  cov[MAXSAT][3]; /* satellite position covariance (m^2) */
    float  vco[MAXSAT][3]; /* satellite velocity covariance (m^2) */
} peph_t;

typedef struct {        /* precise clock type */
    gtime_t time;       /* time (GPST) */
    int index;          /* clock index for multiple files */
    double clk[MAXSAT][1]; /* satellite clock (s) */
    float  std[MAXSAT][1]; /* satellite clock std (s) */
} pclk_t;

typedef struct {        /* SBAS ephemeris type */
    int sat;            /* satellite number */
    gtime_t t0;         /* reference epoch time (GPST) */
    gtime_t tof;        /* time of message frame (GPST) */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    double pos[3];      /* satellite position (m) (ecef) */
    double vel[3];      /* satellite velocity (m/s) (ecef) */
    double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
    double af0,af1;     /* satellite clock-offset/drift (s,s/s) */
} seph_t;

typedef struct {        /* NORAL TLE data type */
    char name [32];     /* common name */
    char alias[32];     /* alias name */
    char satno[16];     /* satellilte catalog number */
    char satclass;      /* classification */
    char desig[16];     /* international designator */
    gtime_t epoch;      /* element set epoch (UTC) */
    double ndot;        /* 1st derivative of mean motion */
    double nddot;       /* 2st derivative of mean motion */
    double bstar;       /* B* drag term */
    int etype;          /* element set type */
    int eleno;          /* element number */
    double inc;         /* orbit inclination (deg) */
    double OMG;         /* right ascension of ascending node (deg) */
    double ecc;         /* eccentricity */
    double omg;         /* argument of perigee (deg) */
    double M;           /* mean anomaly (deg) */
    double n;           /* mean motion (rev/day) */
    int rev;            /* revolution number at epoch */
} tled_t;

typedef struct {        /* NORAD TLE (two line element) type */
    int n,nmax;         /* number/max number of two line element data */
    tled_t *data;       /* NORAD TLE data */
} tle_t;

typedef struct {        /* TEC grid type */
    gtime_t time;       /* epoch time (GPST) */
    int ndata[3];       /* TEC grid data size {nlat,nlon,nhgt} */
    double rb;          /* earth radius (km) */
    double lats[3];     /* latitude start/interval (deg) */
    double lons[3];     /* longitude start/interval (deg) */
    double hgts[3];     /* heights start/interval (km) */
    double *data;       /* TEC grid data (tecu) */
    float *rms;         /* RMS values (tecu) */
} tec_t;

typedef struct {        /* SBAS message type */
    int week,tow;       /* receiption time */
    uint8_t prn,rcv;    /* SBAS satellite PRN,receiver number */
    uint8_t msg[29];    /* SBAS message (226bit) padded by 0 */
} sbsmsg_t;

typedef struct {        /* SBAS messages type */
    int n,nmax;         /* number of SBAS messages/allocated */
    sbsmsg_t *msgs;     /* SBAS messages */
} sbs_t;

typedef struct {        /* SBAS fast correction type */
    gtime_t t0;         /* time of applicability (TOF) */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range-rate correction (RRC) (m/s) */
    double dt;          /* range-rate correction delta-time (s) */
    int iodf;           /* IODF (issue of date fast corr) */
    int16_t udre;       /* UDRE+1 */
    int16_t ai;         /* degradation factor indicator */
} sbsfcorr_t;

typedef struct {        /* SBAS long term satellite error correction type */
    gtime_t t0;         /* correction time */
    int iode;           /* IODE (issue of date ephemeris) */
    double dpos[3];     /* delta position (m) (ecef) */
    double dvel[3];     /* delta velocity (m/s) (ecef) */
    double daf0,daf1;   /* delta clock-offset/drift (s,s/s) */
} sbslcorr_t;

typedef struct {        /* SBAS satellite correction type */
    int sat;            /* satellite number */
    sbsfcorr_t fcorr;   /* fast correction */
    sbslcorr_t lcorr;   /* long term correction */
} sbssatp_t;

typedef struct {        /* SBAS satellite corrections type */
    int iodp;           /* IODP (issue of date mask) */
    int nsat;           /* number of satellites */
    int tlat;           /* system latency (s) */
    sbssatp_t sat[MAXSAT]; /* satellite correction */
} sbssat_t;

typedef struct {        /* SBAS ionospheric correction type */
    gtime_t t0;         /* correction time */
    int16_t lat,lon;    /* latitude/longitude (deg) */
    int16_t give;       /* GIVI+1 */
    float delay;        /* vertical delay estimate (m) */
} sbsigp_t;

typedef struct {        /* IGP band type */
    int16_t x;          /* longitude/latitude (deg) */
    const int16_t *y;   /* latitudes/longitudes (deg) */
    uint8_t bits;       /* IGP mask start bit */
    uint8_t bite;       /* IGP mask end bit */
} sbsigpband_t;

typedef struct {        /* SBAS ionospheric corrections type */
    int iodi;           /* IODI (issue of date ionos corr) */
    int nigp;           /* number of igps */
    sbsigp_t igp[MAXNIGP]; /* ionospheric correction */
} sbsion_t;

typedef struct {        /* DGPS/GNSS correction type */
    gtime_t t0;         /* correction time */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range rate correction (RRC) (m/s) */
    int iod;            /* issue of data (IOD) */
    double udre;        /* UDRE */
} dgps_t;

typedef struct {        /* SSR correction type */
    gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    double udi[6];      /* SSR update interval (s) */
    int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;           /* issue of data */
    int iodcrc;         /* issue of data crc for beidou/sbas */
    int ura;            /* URA indicator */
    int refd;           /* sat ref datum (0:ITRF,1:regional) */
    double deph [3];    /* delta orbit {radial,along,cross} (m) */
    double ddeph[3];    /* dot delta orbit {radial,along,cross} (m/s) */
    double dclk [3];    /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    double hrclk;       /* high-rate clock corection (m) */
    float  cbias[MAXCODE]; /* code biases (m) */
    double pbias[MAXCODE]; /* phase biases (m) */
    float  stdpb[MAXCODE]; /* std-dev of phase biases (m) */
    double yaw_ang,yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
    uint8_t update;     /* update flag (0:no update,1:update) */
} ssr_t;

typedef struct {        /* navigation data type */
    int n,nmax;         /* number of broadcast ephemeris */
    int ng,ngmax;       /* number of glonass ephemeris */
    int ns,nsmax;       /* number of sbas ephemeris */
    int ne,nemax;       /* number of precise ephemeris */
    int nc,ncmax;       /* number of precise clock */
    int na,namax;       /* number of almanac data */
    int nt,ntmax;       /* number of tec grid data */
    eph_t *eph;         /* GPS/QZS/GAL/BDS/IRN ephemeris */
    geph_t *geph;       /* GLONASS ephemeris */
    seph_t *seph;       /* SBAS ephemeris */
    peph_t *peph;       /* precise ephemeris */
    pclk_t *pclk;       /* precise clock */
    alm_t *alm;         /* almanac data */
    tec_t *tec;         /* tec grid data */
    erp_t  erp;         /* earth rotation parameters */
    double utc_gps[8];  /* GPS delta-UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
    double utc_glo[8];  /* GLONASS UTC time parameters {tau_C,tau_GPS} */
    double utc_gal[8];  /* Galileo UTC parameters */
    double utc_qzs[8];  /* QZS UTC parameters */
    double utc_cmp[8];  /* BeiDou UTC parameters */
    double utc_irn[9];  /* IRNSS UTC parameters {A0,A1,Tot,...,dt_LSF,A2} */
    double utc_sbs[4];  /* SBAS UTC parameters */
    double ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int glo_fcn[32];    /* GLONASS FCN + 8 */
    double cbias[MAXSAT][3]; /* satellite DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
    double rbias[MAXRCV][2][3]; /* receiver DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
    pcv_t pcvs[MAXSAT]; /* satellite antenna pcv */
    sbssat_t sbssat;    /* SBAS satellite corrections */
    sbsion_t sbsion[MAXBAND+1]; /* SBAS ionosphere corrections */
    dgps_t dgps[MAXSAT]; /* DGPS corrections */
    ssr_t ssr[MAXSAT];  /* SSR corrections */
} nav_t;

typedef struct {        /* station parameter type */
    char name   [MAXANT]; /* marker name */
    char marker [MAXANT]; /* marker number */
    char antdes [MAXANT]; /* antenna descriptor */
    char antsno [MAXANT]; /* antenna serial number */
    char rectype[MAXANT]; /* receiver type descriptor */
    char recver [MAXANT]; /* receiver firmware version */
    char recsno [MAXANT]; /* receiver serial number */
    int antsetup;       /* antenna setup id */
    int itrf;           /* ITRF realization year */
    int deltype;        /* antenna delta type (0:enu,1:xyz) */
    double pos[3];      /* station position (ecef) (m) */
    double del[3];      /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;         /* antenna height (m) */
    int glo_cp_align;   /* GLONASS code-phase alignment (0:no,1:yes) */
    double glo_cp_bias[4]; /* GLONASS code-phase biases {1C,1P,2C,2P} (m) */
} sta_t;

typedef struct {        /* solution type */
    gtime_t time;       /* time (GPST) */
    double rr[6];       /* position/velocity (m|m/s) */
                        /* {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    float  qr[6];       /* position variance/covariance (m^2) */
                        /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                        /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    float  qv[6];       /* velocity variance/covariance (m^2/s^2) */
    double dtr[6];      /* receiver clock bias to time systems (s) */
    uint8_t type;       /* type (0:xyz-ecef,1:enu-baseline) */
    uint8_t stat;       /* solution status (SOLQ_???) */
    uint8_t ns;         /* number of valid satellites */
    float age;          /* age of differential (s) */
    float ratio;        /* AR ratio factor for valiation */
    float thres;        /* AR ratio threshold for valiation */
} sol_t;

typedef struct {        /* solution buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    int cyclic;         /* cyclic buffer flag */
    int start,end;      /* start/end index */
    gtime_t time;       /* current solution time */
    sol_t *data;        /* solution data */
    double rb[3];       /* reference position {x,y,z} (ecef) (m) */
    uint8_t buff[MAXSOLMSG+1]; /* message buffer */
    int nb;             /* number of byte in message buffer */
} solbuf_t;

typedef struct {        /* solution status type */
    gtime_t time;       /* time (GPST) */
    uint8_t sat;        /* satellite number */
    uint8_t frq;        /* frequency (1:L1,2:L2,...) */
    float az,el;        /* azimuth/elevation angle (rad) */
    float resp;         /* pseudorange residual (m) */
    float resc;         /* carrier-phase residual (m) */
    uint8_t flag;       /* flags: (vsat<<5)+(slip<<3)+fix */
    uint16_t snr;       /* signal strength (*SNR_UNIT dBHz) */
    uint16_t lock;      /* lock counter */
    uint16_t outc;      /* outage counter */
    uint16_t slipc;     /* slip counter */
    uint16_t rejc;      /* reject counter */
} solstat_t;

typedef struct {        /* solution status buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    solstat_t *data;    /* solution status data */
} solstatbuf_t;

typedef struct {        /* RTCM control struct type */
    int staid;          /* station id */
    int stah;           /* station health */
    int seqno;          /* sequence number for rtcm 2 or iods msm */
    int outtype;        /* output message type */
    gtime_t time;       /* message time */
    gtime_t time_s;     /* message start time */
    obs_t obs;          /* observation data (uncorrected) */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    dgps_t *dgps;       /* output of dgps corrections */
    ssr_t ssr[MAXSAT];  /* output of ssr corrections */
    char msg[128];      /* special message */
    char msgtype[256];  /* last message type */
    char msmtype[7][128]; /* msm signal types */
    int obsflag;        /* obs data complete flag (1:ok,0:not complete) */
    int ephsat;         /* input ephemeris satellite number */
    int ephset;         /* input ephemeris set (0-1) */
    double cp[MAXSAT][NFREQ+NEXOBS]; /* carrier-phase measurement */
    uint16_t lock[MAXSAT][NFREQ+NEXOBS]; /* lock time */
    uint16_t loss[MAXSAT][NFREQ+NEXOBS]; /* loss of lock count */
    gtime_t lltime[MAXSAT][NFREQ+NEXOBS]; /* last lock time */
    int nbyte;          /* number of bytes in message buffer */ 
    int nbit;           /* number of bits in word buffer */ 
    int len;            /* message length (bytes) */
    uint8_t buff[1200]; /* message buffer */
    uint32_t word;      /* word buffer for rtcm 2 */
    uint32_t nmsg2[100]; /* message count of RTCM 2 (1-99:1-99,0:other) */
    uint32_t nmsg3[400]; /* message count of RTCM 3 (1-299:1001-1299,300-329:4070-4099,0:ohter) */
    char opt[256];      /* RTCM dependent options */
} rtcm_t;

typedef struct {        /* RINEX control struct type */
    gtime_t time;       /* message time */
    double ver;         /* RINEX version */
    char   type;        /* RINEX file type ('O','N',...) */
    int    sys;         /* navigation system */
    int    tsys;        /* time system */
    char   tobs[8][MAXOBSTYPE][4]; /* rinex obs types */
    obs_t  obs;         /* observation data */
    nav_t  nav;         /* navigation data */
    sta_t  sta;         /* station info */
    int    ephsat;      /* input ephemeris satellite number */
    int    ephset;      /* input ephemeris set (0-1) */
    char   opt[256];    /* rinex dependent options */
} rnxctr_t;

typedef struct {        /* download URL type */
    char type[32];      /* data type */
    char path[1024];    /* URL path */
    char dir [1024];    /* local directory */
    double tint;        /* time interval (s) */
} url_t;

typedef struct {        /* option type */
    const char *name;   /* option name */
    int format;         /* option format (0:int,1:double,2:string,3:enum) */
    void *var;          /* pointer to option variable */
    const char *comment; /* option comment/enum labels/unit */
} opt_t;

typedef struct {        /* SNR mask type */
    int ena[2];         /* enable flag {rover,base} */
    double mask[NFREQ][9]; /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;

typedef struct {        /* processing options type */
    int mode;           /* positioning mode (PMODE_???) */
    int soltype;        /* solution type (0:forward,1:backward,2:combined) */
    int nf;             /* number of frequencies (1:L1,2:L1+L2,3:L1+L2+L5) */
    int navsys;         /* navigation system */
    double elmin;       /* elevation mask angle (rad) */
    snrmask_t snrmask;  /* SNR mask */
    int sateph;         /* satellite ephemeris/clock (EPHOPT_???) */
    int modear;         /* AR mode (0:off,1:continuous,2:instantaneous,3:fix and hold,4:ppp-ar) */
    int glomodear;      /* GLONASS AR mode (0:off,1:on,2:auto cal,3:ext cal) */
    int bdsmodear;      /* BeiDou AR mode (0:off,1:on) */
    int maxout;         /* obs outage count to reset bias */
    int minlock;        /* min lock count to fix ambiguity */
    int minfix;         /* min fix count to hold ambiguity */
    int armaxiter;      /* max iteration to resolve ambiguity */
    int ionoopt;        /* ionosphere option (IONOOPT_???) */
    int tropopt;        /* troposphere option (TROPOPT_???) */
    int dynamics;       /* dynamics model (0:none,1:velociy,2:accel) */
    int tidecorr;       /* earth tide correction (0:off,1:solid,2:solid+otl+pole) */
    int niter;          /* number of filter iteration */
    int codesmooth;     /* code smoothing window size (0:none) */
    int intpref;        /* interpolate reference obs (for post mission) */
    int sbascorr;       /* SBAS correction options */
    int sbassatsel;     /* SBAS satellite selection (0:all) */
    int rovpos;         /* rover position for fixed mode */
    int refpos;         /* base position for relative mode */
                        /* (0:pos in prcopt,  1:average of single pos, */
                        /*  2:read from file, 3:rinex header, 4:rtcm pos) */
    double eratio[NFREQ]; /* code/phase error ratio */
    double err[5];      /* measurement error factor */
                        /* [0]:reserved */
                        /* [1-3]:error factor a/b/c of phase (m) */
                        /* [4]:doppler frequency (hz) */
    double std[3];      /* initial-state std [0]bias,[1]iono [2]trop */
    double prn[6];      /* process-noise std [0]bias,[1]iono [2]trop [3]acch [4]accv [5] pos */
    double sclkstab;    /* satellite clock stability (sec/sec) */
    double thresar[8];  /* AR validation threshold */
    double elmaskar;    /* elevation mask of AR for rising satellite (deg) */
    double elmaskhold;  /* elevation mask to hold ambiguity (deg) */
    double thresslip;   /* slip threshold of geometry-free phase (m) */
    double maxtdiff;    /* max difference of time (sec) */
    double maxinno;     /* reject threshold of innovation (m) */
    double maxgdop;     /* reject threshold of gdop */
    double baseline[2]; /* baseline length constraint {const,sigma} (m) */
    double ru[3];       /* rover position for fixed mode {x,y,z} (ecef) (m) */
    double rb[3];       /* base position for relative mode {x,y,z} (ecef) (m) */
    char anttype[2][MAXANT]; /* antenna types {rover,base} */
    double antdel[2][3]; /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
    pcv_t pcvr[2];      /* receiver antenna parameters {rov,base} */
    uint8_t exsats[MAXSAT]; /* excluded satellites (1:excluded,2:included) */
    int  maxaveep;      /* max averaging epoches */
    int  initrst;       /* initialize by restart */
    int  outsingle;     /* output single by dgps/float/fix/ppp outage */
    char rnxopt[2][256]; /* rinex options {rover,base} */
    int  posopt[6];     /* positioning options */
    int  syncsol;       /* solution sync mode (0:off,1:on) */
    double odisp[2][6*11]; /* ocean tide loading parameters {rov,base} */
    int  freqopt;       /* disable L2-AR */
    char pppopt[256];   /* ppp option */
} prcopt_t;

typedef struct {        /* solution options type */
    int posf;           /* solution format (SOLF_???) */
    int times;          /* time system (TIMES_???) */
    int timef;          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;          /* time digits under decimal point */
    int degf;           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outhead;        /* output header (0:no,1:yes) */
    int outopt;         /* output processing options (0:no,1:yes) */
    int outvel;         /* output velocity options (0:no,1:yes) */
    int datum;          /* datum (0:WGS84,1:Tokyo) */
    int height;         /* height (0:ellipsoidal,1:geodetic) */
    int geoid;          /* geoid model (0:EGM96,1:JGD2000) */
    int solstatic;      /* solution of static mode (0:all,1:single) */
    int sstat;          /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;          /* debug trace level (0:off,1-5:debug) */
    double nmeaintv[2]; /* nmea output interval (s) (<0:no,0:all) */
                        /* nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
    char sep[64];       /* field separator */
    char prog[64];      /* program name */
    double maxsolstd;   /* max std-dev for solution output (m) (0:all) */
} solopt_t;

typedef struct {        /* file options type */
    char satantp[MAXSTRPATH]; /* satellite antenna parameters file */
    char rcvantp[MAXSTRPATH]; /* receiver antenna parameters file */
    char stapos [MAXSTRPATH]; /* station positions file */
    char geoid  [MAXSTRPATH]; /* external geoid data file */
    char iono   [MAXSTRPATH]; /* ionosphere data file */
    char dcb    [MAXSTRPATH]; /* dcb data file */
    char eop    [MAXSTRPATH]; /* eop data file */
    char blq    [MAXSTRPATH]; /* ocean tide loading blq file */
    char tempdir[MAXSTRPATH]; /* ftp/http temporaly directory */
    char geexe  [MAXSTRPATH]; /* google earth exec file */
    char solstat[MAXSTRPATH]; /* solution statistics file */
    char trace  [MAXSTRPATH]; /* debug trace file */
} filopt_t;

typedef struct {        /* RINEX options type */
    gtime_t ts,te;      /* time start/end */
    double tint;        /* time interval (s) */
    double ttol;        /* time tolerance (s) */
    double tunit;       /* time unit for multiple-session (s) */
    int rnxver;         /* RINEX version (x100) */
    int navsys;         /* navigation system */
    int obstype;        /* observation type */
    int freqtype;       /* frequency type */
    char mask[7][64];   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    char staid [32];    /* station id for rinex file name */
    char prog  [32];    /* program */
    char runby [32];    /* run-by */
    char marker[64];    /* marker name */
    char markerno[32];  /* marker number */
    char markertype[32]; /* marker type (ver.3) */
    char name[2][32];   /* observer/agency */
    char rec [3][32];   /* receiver #/type/vers */
    char ant [3][32];   /* antenna #/type */
    double apppos[3];   /* approx position x/y/z */
    double antdel[3];   /* antenna delta h/e/n */
    double glo_cp_bias[4]; /* GLONASS code-phase biases (m) */
    char comment[MAXCOMMENT][64]; /* comments */
    char rcvopt[256];   /* receiver dependent options */
    uint8_t exsats[MAXSAT]; /* excluded satellites */
    int glofcn[32];     /* glonass fcn+8 */
    int outiono;        /* output iono correction */
    int outtime;        /* output time system correction */
    int outleaps;       /* output leap seconds */
    int autopos;        /* auto approx position */
    int phshift;        /* phase shift correction */
    int halfcyc;        /* half cycle correction */
    int sep_nav;        /* separated nav files */
    gtime_t tstart;     /* first obs time */
    gtime_t tend;       /* last obs time */
    gtime_t trtcm;      /* approx log start time for rtcm */
    char tobs[7][MAXOBSTYPE][4]; /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    double shift[7][MAXOBSTYPE]; /* phase shift (cyc) {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    int nobs[7];        /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
} rnxopt_t;

typedef struct {        /* satellite status type */
    uint8_t sys;        /* navigation system */
    uint8_t vs;         /* valid satellite flag single */
    double azel[2];     /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ]; /* residuals of pseudorange (m) */
    double resc[NFREQ]; /* residuals of carrier-phase (m) */
    uint8_t vsat[NFREQ]; /* valid satellite flag */
    uint16_t snr[NFREQ]; /* signal strength (*SNR_UNIT dBHz) */
    uint8_t fix [NFREQ]; /* ambiguity fix flag (1:fix,2:float,3:hold) */
    uint8_t slip[NFREQ]; /* cycle-slip flag */
    uint8_t half[NFREQ]; /* half-cycle valid flag */
    int lock [NFREQ];   /* lock counter of phase */
    uint32_t outc [NFREQ]; /* obs outage counter of phase */
    uint32_t slipc[NFREQ]; /* cycle-slip counter */
    uint32_t rejc [NFREQ]; /* reject counter */
    double gf[NFREQ-1]; /* geometry-free phase (m) */
    double mw[NFREQ-1]; /* MW-LC (m) */
    double phw;         /* phase windup (cycle) */
    gtime_t pt[2][NFREQ]; /* previous carrier-phase time */
    double ph[2][NFREQ]; /* previous carrier-phase observable (cycle) */
} ssat_t;

typedef struct {        /* ambiguity control type */
    gtime_t epoch[4];   /* last epoch */
    int n[4];           /* number of epochs */
    double LC [4];      /* linear combination average */
    double LCv[4];      /* linear combination variance */
    int fixcnt;         /* fix count */
    char flags[MAXSAT]; /* fix flags */
} ambc_t;

typedef struct {        /* RTK control/result type */
    sol_t  sol;         /* RTK solution */
    double rb[6];       /* base position/velocity (ecef) (m|m/s) */
    int nx,na;          /* number of float states/fixed states */
    double tt;          /* time difference between current and previous (s) */
    double *x, *P;      /* float states and their covariance */
    double *xa,*Pa;     /* fixed states and their covariance */
    int nfix;           /* number of continuous fixes of ambiguity */
    ambc_t ambc[MAXSAT]; /* ambibuity control */
    ssat_t ssat[MAXSAT]; /* satellite status */
    int neb;            /* bytes in error message buffer */
    char errbuf[MAXERRMSG]; /* error message buffer */
    prcopt_t opt;       /* processing options */
} rtk_t;

typedef struct {        /* receiver raw data control type */
    gtime_t time;       /* message time */
    gtime_t tobs[MAXSAT][NFREQ+NEXOBS]; /* observation data time */
    obs_t obs;          /* observation data */
    obs_t obuf;         /* observation data buffer */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    int ephsat;         /* update satelle of ephemeris (0:no satellite) */
    int ephset;         /* update set of ephemeris (0-1) */
    sbsmsg_t sbsmsg;    /* SBAS message */
    char msgtype[256];  /* last message type */
    uint8_t subfrm[MAXSAT][380]; /* subframe buffer */
    double lockt[MAXSAT][NFREQ+NEXOBS]; /* lock time (s) */
    double icpp[MAXSAT],off[MAXSAT],icpc; /* carrier params for ss2 */
    double prCA[MAXSAT],dpCA[MAXSAT]; /* L1/CA pseudrange/doppler for javad */
    uint8_t halfc[MAXSAT][NFREQ+NEXOBS]; /* half-cycle add flag */
    char freqn[MAXOBS]; /* frequency number for javad */
    int nbyte;          /* number of bytes in message buffer */ 
    int len;            /* message length (bytes) */
    int iod;            /* issue of data */
    int tod;            /* time of day (ms) */
    int tbase;          /* time base (0:gpst,1:utc(usno),2:glonass,3:utc(su) */
    int flag;           /* general purpose flag */
    int outtype;        /* output message type */
    uint8_t buff[MAXRAWLEN]; /* message buffer */
    char opt[256];      /* receiver dependent options */
    int format;         /* receiver stream format */
    void *rcv_data;     /* receiver dependent data */
} raw_t;

typedef struct {        /* stream type */
    int type;           /* type (STR_???) */
    int mode;           /* mode (STR_MODE_?) */
    int state;          /* state (-1:error,0:close,1:open) */
    uint32_t inb,inr;   /* input bytes/rate */
    uint32_t outb,outr; /* output bytes/rate */
    uint32_t tick_i;    /* input tick tick */
    uint32_t tick_o;    /* output tick */
    uint32_t tact;      /* active tick */
    uint32_t inbt,outbt; /* input/output bytes at tick */
    lock_t lock;        /* lock flag */
    void *port;         /* type dependent port control struct */
    char path[MAXSTRPATH]; /* stream path */
    char msg [MAXSTRMSG];  /* stream message */
} stream_t;

typedef struct {        /* stream converter type */
    int itype,otype;    /* input and output stream type */
    int nmsg;           /* number of output messages */
    int msgs[32];       /* output message types */
    double tint[32];    /* output message intervals (s) */
    uint32_t tick[32];  /* cycle tick of output message */
    int ephsat[32];     /* satellites of output ephemeris */
    int stasel;         /* station info selection (0:remote,1:local) */
    rtcm_t rtcm;        /* rtcm input data buffer */
    raw_t raw;          /* raw  input data buffer */
    rtcm_t out;         /* rtcm output data buffer */
} strconv_t;

typedef struct {        /* stream server type */
    int state;          /* server state (0:stop,1:running) */
    int cycle;          /* server cycle (ms) */
    int buffsize;       /* input/monitor buffer size (bytes) */
    int nmeacycle;      /* NMEA request cycle (ms) (0:no) */
    int relayback;      /* relay back of output streams (0:no) */
    int nstr;           /* number of streams (1 input + (nstr-1) outputs */
    int npb;            /* data length in peek buffer (bytes) */
    char cmds_periodic[16][MAXRCVCMD]; /* periodic commands */
    double nmeapos[3];  /* NMEA request position (ecef) (m) */
    uint8_t *buff;      /* input buffers */
    uint8_t *pbuf;      /* peek buffer */
    uint32_t tick;      /* start tick */
    stream_t stream[16]; /* input/output streams */
    stream_t strlog[16]; /* return log streams */
    strconv_t *conv[16]; /* stream converter */
    thread_t thread;    /* server thread */
    lock_t lock;        /* lock flag */
} strsvr_t;

typedef struct {        /* RTK server type */
    int state;          /* server state (0:stop,1:running) */
    int cycle;          /* processing cycle (ms) */
    int nmeacycle;      /* NMEA request cycle (ms) (0:no req) */
    int nmeareq;        /* NMEA request (0:no,1:nmeapos,2:single sol) */
    double nmeapos[3];  /* NMEA request position (ecef) (m) */
    int buffsize;       /* input buffer size (bytes) */
    int format[3];      /* input format {rov,base,corr} */
    solopt_t solopt[2]; /* output solution options {sol1,sol2} */
    int navsel;         /* ephemeris select (0:all,1:rover,2:base,3:corr) */
    int nsbs;           /* number of sbas message */
    int nsol;           /* number of solution buffer */
    rtk_t rtk;          /* RTK control/result struct */
    int nb [3];         /* bytes in input buffers {rov,base} */
    int nsb[2];         /* bytes in soulution buffers */
    int npb[3];         /* bytes in input peek buffers */
    uint8_t *buff[3];   /* input buffers {rov,base,corr} */
    uint8_t *sbuf[2];   /* output buffers {sol1,sol2} */
    uint8_t *pbuf[3];   /* peek buffers {rov,base,corr} */
    sol_t solbuf[MAXSOLBUF]; /* solution buffer */
    uint32_t nmsg[3][10]; /* input message counts */
    raw_t  raw [3];     /* receiver raw control {rov,base,corr} */
    rtcm_t rtcm[3];     /* RTCM control {rov,base,corr} */
    gtime_t ftime[3];   /* download time {rov,base,corr} */
    char files[3][MAXSTRPATH]; /* download paths {rov,base,corr} */
    obs_t obs[3][MAXOBSBUF]; /* observation data {rov,base,corr} */
    nav_t nav;          /* navigation data */
    sbsmsg_t sbsmsg[MAXSBSMSG]; /* SBAS message buffer */
    stream_t stream[8]; /* streams {rov,base,corr,sol1,sol2,logr,logb,logc} */
    stream_t *moni;     /* monitor stream */
    uint32_t tick;      /* start tick */
    thread_t thread;    /* server thread */
    int cputime;        /* CPU time (ms) for a processing cycle */
    int prcout;         /* missing observation data count */
    int nave;           /* number of averaging base pos */
    double rb_ave[3];   /* averaging base pos */
    char cmds_periodic[3][MAXRCVCMD]; /* periodic commands */
    char cmd_reset[MAXRCVCMD]; /* reset command */
    double bl_reset;    /* baseline length to reset (km) */
    lock_t lock;        /* lock flag */
} rtksvr_t;

typedef struct {        /* GIS data point type */
    double pos[3];      /* point data {lat,lon,height} (rad,m) */
} gis_pnt_t;

typedef struct {        /* GIS data polyline type */
    int npnt;           /* number of points */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
    double *pos;        /* position data (3 x npnt) */
} gis_poly_t;

typedef struct {        /* GIS data polygon type */
    int npnt;           /* number of points */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
    double *pos;        /* position data (3 x npnt) */
} gis_polygon_t;

typedef struct gisd_tag { /* GIS data list type */
    int type;           /* data type (1:point,2:polyline,3:polygon) */
    void *data;         /* data body */
    struct gisd_tag *next; /* pointer to next */
} gisd_t;

typedef struct {        /* GIS type */
    char name[MAXGISLAYER][256]; /* name */
    int flag[MAXGISLAYER];     /* flag */
    gisd_t *data[MAXGISLAYER]; /* gis data list */
    double bound[4];    /* boundary {lat0,lat1,lon0,lon1} */
} gis_t;

typedef void fatalfunc_t(const char *); /* fatal callback function type */

/* global variables ----------------------------------------------------------*/
extern const double chisqr[];        /* chi-sqr(n) table (alpha=0.001) */
extern const prcopt_t prcopt_default; /* default positioning options */
extern const solopt_t solopt_default; /* default solution output options */
extern const sbsigpband_t igpband1[9][8]; /* SBAS IGP band 0-8 */
extern const sbsigpband_t igpband2[2][5]; /* SBAS IGP band 9-10 */
extern const char *formatstrs[];     /* stream format strings */
extern opt_t sysopts[];              /* system options table */

/* satellites, systems, codes functions --------------------------------------*/
EXPORT int  satno   (int sys, int prn);
EXPORT int  satsys  (int sat, int *prn);
EXPORT int  satid2no(const char *id);
EXPORT void satno2id(int sat, char *id);
EXPORT uint8_t obs2code(const char *obs);
EXPORT char *code2obs(uint8_t code);
EXPORT double code2freq(int sys, uint8_t code, int fcn);
EXPORT double sat2freq(int sat, uint8_t code, const nav_t *nav);
EXPORT int  code2idx(int sys, uint8_t code);
EXPORT int  satexclude(int sat, double var, int svh, const prcopt_t *opt);
EXPORT int  testsnr(int base, int freq, double el, double snr,
                    const snrmask_t *mask);
EXPORT void setcodepri(int sys, int idx, const char *pri);
EXPORT int  getcodepri(int sys, uint8_t code, const char *opt);

/* matrix and vector functions -----------------------------------------------*/
EXPORT double *mat  (int n, int m);
EXPORT int    *imat (int n, int m);
EXPORT double *zeros(int n, int m);
EXPORT double *eye  (int n);
EXPORT double dot (const double *a, const double *b, int n);
EXPORT double norm(const double *a, int n);
EXPORT void cross3(const double *a, const double *b, double *c);
EXPORT int  normv3(const double *a, double *b);
EXPORT void matcpy(double *A, const double *B, int n, int m);
EXPORT void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
EXPORT int  matinv(double *A, int n);
EXPORT int  solve (const char *tr, const double *A, const double *Y, int n,
                   int m, double *X);
EXPORT int  lsq   (const double *A, const double *y, int n, int m, double *x,
                   double *Q);
EXPORT int  filter(double *x, double *P, const double *H, const double *v,
                   const double *R, int n, int m);
EXPORT int  smoother(const double *xf, const double *Qf, const double *xb,
                     const double *Qb, int n, double *xs, double *Qs);
EXPORT void matprint (const double *A, int n, int m, int p, int q);
EXPORT void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);

EXPORT void add_fatal(fatalfunc_t *func);

/* time and string functions -------------------------------------------------*/
EXPORT double  str2num(const char *s, int i, int n);
EXPORT int     str2time(const char *s, int i, int n, gtime_t *t);
EXPORT void    time2str(gtime_t t, char *str, int n);
EXPORT gtime_t epoch2time(const double *ep);
EXPORT void    time2epoch(gtime_t t, double *ep);
EXPORT gtime_t gpst2time(int week, double sec);
EXPORT double  time2gpst(gtime_t t, int *week);
EXPORT gtime_t gst2time(int week, double sec);
EXPORT double  time2gst(gtime_t t, int *week);
EXPORT gtime_t bdt2time(int week, double sec);
EXPORT double  time2bdt(gtime_t t, int *week);
EXPORT char    *time_str(gtime_t t, int n);

EXPORT gtime_t timeadd  (gtime_t t, double sec);
EXPORT double  timediff (gtime_t t1, gtime_t t2);
EXPORT gtime_t gpst2utc (gtime_t t);
EXPORT gtime_t utc2gpst (gtime_t t);
EXPORT gtime_t gpst2bdt (gtime_t t);
EXPORT gtime_t bdt2gpst (gtime_t t);
EXPORT gtime_t timeget  (void);
EXPORT void    timeset  (gtime_t t);
EXPORT void    timereset(void);
EXPORT double  time2doy (gtime_t t);
EXPORT double  utc2gmst (gtime_t t, double ut1_utc);
EXPORT int read_leaps(const char *file);

EXPORT int adjgpsweek(int week);
EXPORT uint32_t tickget(void);
EXPORT void sleepms(int ms);

EXPORT int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
                   const char *base);
EXPORT int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts,
                    gtime_t te, const char *rov, const char *base);

/* coordinates transformation ------------------------------------------------*/
EXPORT void ecef2pos(const double *r, double *pos);
EXPORT void pos2ecef(const double *pos, double *r);
EXPORT void ecef2enu(const double *pos, const double *r, double *e);
EXPORT void enu2ecef(const double *pos, const double *e, double *r);
EXPORT void covenu  (const double *pos, const double *P, double *Q);
EXPORT void covecef (const double *pos, const double *Q, double *P);
EXPORT void xyz2enu (const double *pos, double *E);
EXPORT void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
EXPORT void deg2dms (double deg, double *dms, int ndec);
EXPORT double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
EXPORT void readpos(const char *file, const char *rcv, double *pos);
EXPORT int  sortobs(obs_t *obs);
EXPORT void uniqnav(nav_t *nav);
EXPORT int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
EXPORT int  readnav(const char *file, nav_t *nav);
EXPORT int  savenav(const char *file, const nav_t *nav);
EXPORT void freeobs(obs_t *obs);
EXPORT void freenav(nav_t *nav, int opt);
EXPORT int  readblq(const char *file, const char *sta, double *odisp);
EXPORT int  readerp(const char *file, erp_t *erp);
EXPORT int  geterp (const erp_t *erp, gtime_t time, double *val);

/* debug trace functions -----------------------------------------------------*/
EXPORT void traceopen(const char *file);
EXPORT void traceclose(void);
EXPORT void tracelevel(int level);
EXPORT void trace    (int level, const char *format, ...);
EXPORT void tracet   (int level, const char *format, ...);
EXPORT void tracemat (int level, const double *A, int n, int m, int p, int q);
EXPORT void traceobs (int level, const obsd_t *obs, int n);
EXPORT void tracenav (int level, const nav_t *nav);
EXPORT void tracegnav(int level, const nav_t *nav);
EXPORT void tracehnav(int level, const nav_t *nav);
EXPORT void tracepeph(int level, const nav_t *nav);
EXPORT void tracepclk(int level, const nav_t *nav);
EXPORT void traceb   (int level, const uint8_t *p, int n);

/* platform dependent functions ----------------------------------------------*/
EXPORT int execcmd(const char *cmd);
EXPORT int expath (const char *path, char *paths[], int nmax);
EXPORT void createdir(const char *path);

/* positioning models --------------------------------------------------------*/
EXPORT double satazel(const double *pos, const double *e, double *azel);
EXPORT double geodist(const double *rs, const double *rr, double *e);
EXPORT void dops(int ns, const double *azel, double elmin, double *dop);

/* atmosphere models ---------------------------------------------------------*/
EXPORT double ionmodel(gtime_t t, const double *ion, const double *pos,
                       const double *azel);
EXPORT double ionmapf(const double *pos, const double *azel);
EXPORT double ionppp(const double *pos, const double *azel, double re,
                     double hion, double *pppos);
EXPORT double tropmodel(gtime_t time, const double *pos, const double *azel,
                        double humi);
EXPORT double tropmapf(gtime_t time, const double *pos, const double *azel,
                       double *mapfw);
EXPORT int iontec(gtime_t time, const nav_t *nav, const double *pos,
                  const double *azel, int opt, double *delay, double *var);
EXPORT void readtec(const char *file, nav_t *nav, int opt);
EXPORT int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var);
EXPORT int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var);

/* antenna models ------------------------------------------------------------*/
EXPORT int  readpcv(const char *file, pcvs_t *pcvs);
EXPORT pcv_t *searchpcv(int sat, const char *type, gtime_t time,
                        const pcvs_t *pcvs);
EXPORT void antmodel(const pcv_t *pcv, const double *del, const double *azel,
                     int opt, double *dant);
EXPORT void antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
EXPORT void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
                       double *rmoon, double *gmst);
EXPORT void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr);

/* geiod models --------------------------------------------------------------*/
EXPORT int opengeoid(int model, const char *file);
EXPORT void closegeoid(void);
EXPORT double geoidh(const double *pos);

/* datum transformation ------------------------------------------------------*/
EXPORT int loaddatump(const char *file);
EXPORT int tokyo2jgd(double *pos);
EXPORT int jgd2tokyo(double *pos);

/* rinex functions -----------------------------------------------------------*/
EXPORT int readrnx (const char *file, int rcv, const char *opt, obs_t *obs,
                    nav_t *nav, sta_t *sta);
EXPORT int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
                    double tint, const char *opt, obs_t *obs, nav_t *nav,
                    sta_t *sta);
EXPORT int readrnxc(const char *file, nav_t *nav);
EXPORT int outrnxobsh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxobsb(FILE *fp, const rnxopt_t *opt, const obsd_t *obs, int n,
                      int epflag);
EXPORT int outrnxnavh (FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxgnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxhnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxlnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxqnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxcnavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxinavh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
EXPORT int outrnxnavb (FILE *fp, const rnxopt_t *opt, const eph_t *eph);
EXPORT int outrnxgnavb(FILE *fp, const rnxopt_t *opt, const geph_t *geph);
EXPORT int outrnxhnavb(FILE *fp, const rnxopt_t *opt, const seph_t *seph);
EXPORT int rtk_uncompress(const char *file, char *uncfile);
EXPORT int convrnx(int format, rnxopt_t *opt, const char *file, char **ofile);
EXPORT int  init_rnxctr (rnxctr_t *rnx);
EXPORT void free_rnxctr (rnxctr_t *rnx);
EXPORT int  open_rnxctr (rnxctr_t *rnx, FILE *fp);
EXPORT int  input_rnxctr(rnxctr_t *rnx, FILE *fp);

/* ephemeris and clock functions ---------------------------------------------*/
EXPORT double eph2clk (gtime_t time, const eph_t  *eph);
EXPORT double geph2clk(gtime_t time, const geph_t *geph);
EXPORT double seph2clk(gtime_t time, const seph_t *seph);
EXPORT void eph2pos (gtime_t time, const eph_t  *eph,  double *rs, double *dts,
                     double *var);
EXPORT void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var);
EXPORT void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var);
EXPORT int  peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                     double *rs, double *dts, double *var);
EXPORT void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant);
EXPORT int  satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                   const nav_t *nav, double *rs, double *dts, double *var,
                   int *svh);
EXPORT void satposs(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                    int sateph, double *rs, double *dts, double *var, int *svh);
EXPORT void setseleph(int sys, int sel);
EXPORT int  getseleph(int sys);
EXPORT void readsp3(const char *file, nav_t *nav, int opt);
EXPORT int  readsap(const char *file, gtime_t time, nav_t *nav);
EXPORT int  readdcb(const char *file, nav_t *nav, const sta_t *sta);
EXPORT int  readfcb(const char *file, nav_t *nav);
EXPORT void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts);

EXPORT int tle_read(const char *file, tle_t *tle);
EXPORT int tle_name_read(const char *file, tle_t *tle);
EXPORT int tle_pos(gtime_t time, const char *name, const char *satno,
                   const char *desig, const tle_t *tle, const erp_t *erp,
                   double *rs);

/* receiver raw data functions -----------------------------------------------*/
EXPORT uint32_t getbitu(const uint8_t *buff, int pos, int len);
EXPORT int32_t  getbits(const uint8_t *buff, int pos, int len);
EXPORT void setbitu(uint8_t *buff, int pos, int len, uint32_t data);
EXPORT void setbits(uint8_t *buff, int pos, int len, int32_t  data);
EXPORT uint32_t rtk_crc32 (const uint8_t *buff, int len);
EXPORT uint32_t rtk_crc24q(const uint8_t *buff, int len);
EXPORT uint16_t rtk_crc16 (const uint8_t *buff, int len);
EXPORT int decode_word (uint32_t word, uint8_t *data);
EXPORT int decode_frame(const uint8_t *buff, eph_t *eph, alm_t *alm,
                        double *ion, double *utc);
EXPORT int test_glostr(const uint8_t *buff);
EXPORT int decode_glostr(const uint8_t *buff, geph_t *geph, double *utc);
EXPORT int decode_bds_d1(const uint8_t *buff, eph_t *eph, double *ion,
                         double *utc);
EXPORT int decode_bds_d2(const uint8_t *buff, eph_t *eph, double *utc);
EXPORT int decode_gal_inav(const uint8_t *buff, eph_t *eph, double *ion,
                           double *utc);
EXPORT int decode_gal_fnav(const uint8_t *buff, eph_t *eph, double *ion,
                           double *utc);
EXPORT int decode_irn_nav(const uint8_t *buff, eph_t *eph, double *ion,
                          double *utc);

EXPORT int init_raw   (raw_t *raw, int format);
EXPORT void free_raw  (raw_t *raw);
EXPORT int input_raw  (raw_t *raw, int format, uint8_t data);
EXPORT int input_rawf (raw_t *raw, int format, FILE *fp);

EXPORT int init_rt17  (raw_t *raw);
EXPORT int init_cmr   (raw_t *raw);
EXPORT void free_rt17 (raw_t *raw);
EXPORT void free_cmr  (raw_t *raw);
EXPORT int update_cmr (raw_t *raw, rtksvr_t *svr, obs_t *obs);

EXPORT int input_oem4  (raw_t *raw, uint8_t data);
EXPORT int input_oem3  (raw_t *raw, uint8_t data);
EXPORT int input_ubx   (raw_t *raw, uint8_t data);
EXPORT int input_ss2   (raw_t *raw, uint8_t data);
EXPORT int input_cres  (raw_t *raw, uint8_t data);
EXPORT int input_stq   (raw_t *raw, uint8_t data);
EXPORT int input_javad (raw_t *raw, uint8_t data);
EXPORT int input_nvs   (raw_t *raw, uint8_t data);
EXPORT int input_bnx   (raw_t *raw, uint8_t data);
EXPORT int input_rt17  (raw_t *raw, uint8_t data);
EXPORT int input_sbf   (raw_t *raw, uint8_t data);
EXPORT int input_oem4f (raw_t *raw, FILE *fp);
EXPORT int input_oem3f (raw_t *raw, FILE *fp);
EXPORT int input_ubxf  (raw_t *raw, FILE *fp);
EXPORT int input_ss2f  (raw_t *raw, FILE *fp);
EXPORT int input_cresf (raw_t *raw, FILE *fp);
EXPORT int input_stqf  (raw_t *raw, FILE *fp);
EXPORT int input_javadf(raw_t *raw, FILE *fp);
EXPORT int input_nvsf  (raw_t *raw, FILE *fp);
EXPORT int input_bnxf  (raw_t *raw, FILE *fp);
EXPORT int input_rt17f (raw_t *raw, FILE *fp);
EXPORT int input_sbff  (raw_t *raw, FILE *fp);

EXPORT int gen_ubx (const char *msg, uint8_t *buff);
EXPORT int gen_stq (const char *msg, uint8_t *buff);
EXPORT int gen_nvs (const char *msg, uint8_t *buff);

/* rtcm functions ------------------------------------------------------------*/
EXPORT int init_rtcm   (rtcm_t *rtcm);
EXPORT void free_rtcm  (rtcm_t *rtcm);
EXPORT int input_rtcm2 (rtcm_t *rtcm, uint8_t data);
EXPORT int input_rtcm3 (rtcm_t *rtcm, uint8_t data);
EXPORT int input_rtcm2f(rtcm_t *rtcm, FILE *fp);
EXPORT int input_rtcm3f(rtcm_t *rtcm, FILE *fp);
EXPORT int gen_rtcm2   (rtcm_t *rtcm, int type, int sync);
EXPORT int gen_rtcm3   (rtcm_t *rtcm, int type, int subtype, int sync);

/* solution functions --------------------------------------------------------*/
EXPORT void initsolbuf(solbuf_t *solbuf, int cyclic, int nmax);
EXPORT void freesolbuf(solbuf_t *solbuf);
EXPORT void freesolstatbuf(solstatbuf_t *solstatbuf);
EXPORT sol_t *getsol(solbuf_t *solbuf, int index);
EXPORT int addsol(solbuf_t *solbuf, const sol_t *sol);
EXPORT int readsol (char *files[], int nfile, solbuf_t *sol);
EXPORT int readsolt(char *files[], int nfile, gtime_t ts, gtime_t te,
                    double tint, int qflag, solbuf_t *sol);
EXPORT int readsolstat(char *files[], int nfile, solstatbuf_t *statbuf);
EXPORT int readsolstatt(char *files[], int nfile, gtime_t ts, gtime_t te,
                        double tint, solstatbuf_t *statbuf);
EXPORT int inputsol(uint8_t data, gtime_t ts, gtime_t te, double tint,
                    int qflag, const solopt_t *opt, solbuf_t *solbuf);

EXPORT int outprcopts(uint8_t *buff, const prcopt_t *opt);
EXPORT int outsolheads(uint8_t *buff, const solopt_t *opt);
EXPORT int outsols  (uint8_t *buff, const sol_t *sol, const double *rb,
                     const solopt_t *opt);
EXPORT int outsolexs(uint8_t *buff, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt);
EXPORT void outprcopt(FILE *fp, const prcopt_t *opt);
EXPORT void outsolhead(FILE *fp, const solopt_t *opt);
EXPORT void outsol  (FILE *fp, const sol_t *sol, const double *rb,
                     const solopt_t *opt);
EXPORT void outsolex(FILE *fp, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt);
EXPORT int outnmea_rmc(uint8_t *buff, const sol_t *sol);
EXPORT int outnmea_gga(uint8_t *buff, const sol_t *sol);
EXPORT int outnmea_gsa(uint8_t *buff, const sol_t *sol,
                       const ssat_t *ssat);
EXPORT int outnmea_gsv(uint8_t *buff, const sol_t *sol,
                       const ssat_t *ssat);

/* google earth kml converter ------------------------------------------------*/
EXPORT int convkml(const char *infile, const char *outfile, gtime_t ts,
                   gtime_t te, double tint, int qflg, double *offset,
                   int tcolor, int pcolor, int outalt, int outtime);

/* gpx converter -------------------------------------------------------------*/
EXPORT int convgpx(const char *infile, const char *outfile, gtime_t ts,
                   gtime_t te, double tint, int qflg, double *offset,
                   int outtrk, int outpnt, int outalt, int outtime);

/* sbas functions ------------------------------------------------------------*/
EXPORT int  sbsreadmsg (const char *file, int sel, sbs_t *sbs);
EXPORT int  sbsreadmsgt(const char *file, int sel, gtime_t ts, gtime_t te,
                        sbs_t *sbs);
EXPORT void sbsoutmsg(FILE *fp, sbsmsg_t *sbsmsg);
EXPORT int  sbsdecodemsg(gtime_t time, int prn, const uint32_t *words,
                         sbsmsg_t *sbsmsg);
EXPORT int sbsupdatecorr(const sbsmsg_t *msg, nav_t *nav);
EXPORT int sbssatcorr(gtime_t time, int sat, const nav_t *nav, double *rs,
                      double *dts, double *var);
EXPORT int sbsioncorr(gtime_t time, const nav_t *nav, const double *pos,
                      const double *azel, double *delay, double *var);
EXPORT double sbstropcorr(gtime_t time, const double *pos, const double *azel,
                          double *var);

/* options functions ---------------------------------------------------------*/
EXPORT opt_t *searchopt(const char *name, const opt_t *opts);
EXPORT int str2opt(opt_t *opt, const char *str);
EXPORT int opt2str(const opt_t *opt, char *str);
EXPORT int opt2buf(const opt_t *opt, char *buff);
EXPORT int loadopts(const char *file, opt_t *opts);
EXPORT int saveopts(const char *file, const char *mode, const char *comment,
                    const opt_t *opts);
EXPORT void resetsysopts(void);
EXPORT void getsysopts(prcopt_t *popt, solopt_t *sopt, filopt_t *fopt);
EXPORT void setsysopts(const prcopt_t *popt, const solopt_t *sopt,
                       const filopt_t *fopt);

/* stream data input and output functions ------------------------------------*/
EXPORT void strinitcom(void);
EXPORT void strinit  (stream_t *stream);
EXPORT void strlock  (stream_t *stream);
EXPORT void strunlock(stream_t *stream);
EXPORT int  stropen  (stream_t *stream, int type, int mode, const char *path);
EXPORT void strclose (stream_t *stream);
EXPORT int  strread  (stream_t *stream, uint8_t *buff, int n);
EXPORT int  strwrite (stream_t *stream, uint8_t *buff, int n);
EXPORT void strsync  (stream_t *stream1, stream_t *stream2);
EXPORT int  strstat  (stream_t *stream, char *msg);
EXPORT int  strstatx (stream_t *stream, char *msg);
EXPORT void strsum   (stream_t *stream, int *inb, int *inr, int *outb, int *outr);
EXPORT void strsetopt(const int *opt);
EXPORT gtime_t strgettime(stream_t *stream);
EXPORT void strsendnmea(stream_t *stream, const sol_t *sol);
EXPORT void strsendcmd(stream_t *stream, const char *cmd);
EXPORT void strsettimeout(stream_t *stream, int toinact, int tirecon);
EXPORT void strsetdir(const char *dir);
EXPORT void strsetproxy(const char *addr);

/* integer ambiguity resolution ----------------------------------------------*/
EXPORT int lambda(int n, int m, const double *a, const double *Q, double *F,
                  double *s);
EXPORT int lambda_reduction(int n, const double *Q, double *Z);
EXPORT int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s);

/* standard positioning ------------------------------------------------------*/
EXPORT int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel,
                  ssat_t *ssat, char *msg);

/* precise positioning -------------------------------------------------------*/
EXPORT void rtkinit(rtk_t *rtk, const prcopt_t *opt);
EXPORT void rtkfree(rtk_t *rtk);
EXPORT int  rtkpos (rtk_t *rtk, const obsd_t *obs, int nobs, const nav_t *nav);
EXPORT int  rtkopenstat(const char *file, int level);
EXPORT void rtkclosestat(void);
EXPORT int  rtkoutstat(rtk_t *rtk, char *buff);

/* precise point positioning -------------------------------------------------*/
EXPORT void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav);
EXPORT int pppnx(const prcopt_t *opt);
EXPORT int pppoutstat(rtk_t *rtk, char *buff);

EXPORT int ppp_ar(rtk_t *rtk, const obsd_t *obs, int n, int *exc,
                  const nav_t *nav, const double *azel, double *x, double *P);

/* post-processing positioning -----------------------------------------------*/
EXPORT int postpos(gtime_t ts, gtime_t te, double ti, double tu,
                   const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, char **infile, int n, char *outfile,
                   const char *rov, const char *base);

/* stream server functions ---------------------------------------------------*/
EXPORT void strsvrinit (strsvr_t *svr, int nout);
EXPORT int  strsvrstart(strsvr_t *svr, int *opts, int *strs, char **paths,
                        char **logs, strconv_t **conv, char **cmds,
                        char **cmds_priodic, const double *nmeapos);
EXPORT void strsvrstop (strsvr_t *svr, char **cmds);
EXPORT void strsvrstat (strsvr_t *svr, int *stat, int *log_stat, int *byte,
                        int *bps, char *msg);
EXPORT strconv_t *strconvnew(int itype, int otype, const char *msgs, int staid,
                             int stasel, const char *opt);
EXPORT void strconvfree(strconv_t *conv);

/* rtk server functions ------------------------------------------------------*/
EXPORT int  rtksvrinit  (rtksvr_t *svr);
EXPORT void rtksvrfree  (rtksvr_t *svr);
EXPORT int  rtksvrstart (rtksvr_t *svr, int cycle, int buffsize, int *strs,
                         char **paths, int *formats, int navsel, char **cmds,
                         char **cmds_periodic, char **rcvopts, int nmeacycle,
                         int nmeareq, const double *nmeapos, prcopt_t *prcopt,
                         solopt_t *solopt, stream_t *moni, char *errmsg);
EXPORT void rtksvrstop  (rtksvr_t *svr, char **cmds);
EXPORT int  rtksvropenstr(rtksvr_t *svr, int index, int str, const char *path,
                          const solopt_t *solopt);
EXPORT void rtksvrclosestr(rtksvr_t *svr, int index);
EXPORT void rtksvrlock  (rtksvr_t *svr);
EXPORT void rtksvrunlock(rtksvr_t *svr);
EXPORT int  rtksvrostat (rtksvr_t *svr, int type, gtime_t *time, int *sat,
                         double *az, double *el, int **snr, int *vsat);
EXPORT void rtksvrsstat (rtksvr_t *svr, int *sstat, char *msg);
EXPORT int  rtksvrmark(rtksvr_t *svr, const char *name, const char *comment);

/* downloader functions ------------------------------------------------------*/
EXPORT int dl_readurls(const char *file, char **types, int ntype, url_t *urls,
                       int nmax);
EXPORT int dl_readstas(const char *file, char **stas, int nmax);
EXPORT int dl_exec(gtime_t ts, gtime_t te, double ti, int seqnos, int seqnoe,
                   const url_t *urls, int nurl, char **stas, int nsta,
                   const char *dir, const char *usr, const char *pwd,
                   const char *proxy, int opts, char *msg, FILE *fp);
EXPORT void dl_test(gtime_t ts, gtime_t te, double ti, const url_t *urls,
                    int nurl, char **stas, int nsta, const char *dir,
                    int ncol, int datefmt, FILE *fp);

/* GIS data functions --------------------------------------------------------*/
EXPORT int gis_read(const char *file, gis_t *gis, int layer);
EXPORT void gis_free(gis_t *gis);

/* application defined functions ---------------------------------------------*/
extern int showmsg(const char *format,...);
extern void settspan(gtime_t ts, gtime_t te);
extern void settime(gtime_t time);

#ifdef __cplusplus
}
#endif
#endif /* RTKLIB_H */


/* ===== Embedded rtkcmn.c ===== */
#pragma push_macro("LAPACK")
#undef LAPACK
#pragma push_macro("MAX_VAR_EPH")
#undef MAX_VAR_EPH
#pragma push_macro("POLYCRC24Q")
#undef POLYCRC24Q
#pragma push_macro("POLYCRC32")
#undef POLYCRC32
#pragma push_macro("Rx")
#undef Rx
#pragma push_macro("Ry")
#undef Ry
#pragma push_macro("Rz")
#undef Rz
#pragma push_macro("SQR")
#undef SQR
#pragma push_macro("SQRT")
#undef SQRT
#pragma push_macro("_POSIX_C_SOURCE")
#undef _POSIX_C_SOURCE
#pragma push_macro("dgemm_")
#undef dgemm_
#pragma push_macro("dgetrf_")
#undef dgetrf_
#pragma push_macro("dgetri_")
#undef dgetri_
#pragma push_macro("dgetrs_")
#undef dgetrs_

/*------------------------------------------------------------------------------
* rtkcmn.c : rtklib common functions
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* options : -DLAPACK   use LAPACK/BLAS
*           -DMKL      use Intel MKL
*           -DTRACE    enable debug trace
*           -DWIN32    use WIN32 API
*           -DNOCALLOC no use calloc for zero matrix
*           -DIERS_MODEL use GMF instead of NMF
*           -DDLL      built for shared library
*           -DCPUTIME_IN_GPST cputime operated in gpst
*
* references :
*     [1] IS-GPS-200D, Navstar GPS Space Segment/Navigation User Interfaces,
*         7 March, 2006
*     [2] RTCA/DO-229C, Minimum operational performance standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         November 28, 2001
*     [3] M.Rothacher, R.Schmid, ANTEX: The Antenna Exchange Format Version 1.4,
*         15 September, 2010
*     [4] A.Gelb ed., Applied Optimal Estimation, The M.I.T Press, 1974
*     [5] A.E.Niell, Global mapping functions for the atmosphere delay at radio
*         wavelengths, Journal of geophysical research, 1996
*     [6] W.Gurtner and L.Estey, RINEX The Receiver Independent Exchange Format
*         Version 3.00, November 28, 2007
*     [7] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*         May 2009
*     [8] China Satellite Navigation Office, BeiDou navigation satellite system
*         signal in space interface control document, open service signal B1I
*         (version 1.0), Dec 2012
*     [9] J.Boehm, A.Niell, P.Tregoning and H.Shuh, Global Mapping Function
*         (GMF): A new empirical mapping function base on numerical weather
*         model data, Geophysical Research Letters, 33, L07304, 2006
*     [10] GLONASS/GPS/Galileo/Compass/SBAS NV08C receiver series BINR interface
*         protocol specification ver.1.3, August, 2012
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/12 1.0 new
*           2007/03/06 1.1 input initial rover pos of pntpos()
*                          update only effective states of filter()
*                          fix bug of atan2() domain error
*           2007/04/11 1.2 add function antmodel()
*                          add gdop mask for pntpos()
*                          change constant MAXDTOE value
*           2007/05/25 1.3 add function execcmd(),expandpath()
*           2008/06/21 1.4 add funciton sortobs(),uniqeph(),screent()
*                          replace geodist() by sagnac correction way
*           2008/10/29 1.5 fix bug of ionospheric mapping function
*                          fix bug of seasonal variation term of tropmapf
*           2008/12/27 1.6 add function tickget(), sleepms(), tracenav(),
*                          xyz2enu(), satposv(), pntvel(), covecef()
*           2009/03/12 1.7 fix bug on error-stop when localtime() returns NULL
*           2009/03/13 1.8 fix bug on time adjustment for summer time
*           2009/04/10 1.9 add function adjgpsweek(),getbits(),getbitu()
*                          add function geph2pos()
*           2009/06/08 1.10 add function seph2pos()
*           2009/11/28 1.11 change function pntpos()
*                           add function tracegnav(),tracepeph()
*           2009/12/22 1.12 change default parameter of ionos std
*                           valid under second for timeget()
*           2010/07/28 1.13 fix bug in tropmapf()
*                           added api:
*                               obs2code(),code2obs(),cross3(),normv3(),
*                               gst2time(),time2gst(),time_str(),timeset(),
*                               deg2dms(),dms2deg(),searchpcv(),antmodel_s(),
*                               tracehnav(),tracepclk(),reppath(),reppaths(),
*                               createdir()
*                           changed api:
*                               readpcv(),
*                           deleted api:
*                               uniqeph()
*           2010/08/20 1.14 omit to include mkl header files
*                           fix bug on chi-sqr(n) table
*           2010/12/11 1.15 added api:
*                               freeobs(),freenav(),ionppp()
*           2011/05/28 1.16 fix bug on half-hour offset by time2epoch()
*                           added api:
*                               uniqnav()
*           2012/06/09 1.17 add a leap second after 2012-6-30
*           2012/07/15 1.18 add api setbits(),setbitu(),utc2gmst()
*                           fix bug on interpolation of antenna pcv
*                           fix bug on str2num() for string with over 256 char
*                           add api readblq(),satexclude(),setcodepri(),
*                           getcodepri()
*                           change api obs2code(),code2obs(),antmodel()
*           2012/12/25 1.19 fix bug on satwavelen(),code2obs(),obs2code()
*                           add api testsnr()
*           2013/01/04 1.20 add api gpst2bdt(),bdt2gpst(),bdt2time(),time2bdt()
*                           readblq(),readerp(),geterp(),crc16()
*                           change api eci2ecef(),sunmoonpos()
*           2013/03/26 1.21 tickget() uses clock_gettime() for linux
*           2013/05/08 1.22 fix bug on nutation coefficients for ast_args()
*           2013/06/02 1.23 add #ifdef for undefined CLOCK_MONOTONIC_RAW
*           2013/09/01 1.24 fix bug on interpolation of satellite antenna pcv
*           2013/09/06 1.25 fix bug on extrapolation of erp
*           2014/04/27 1.26 add SYS_LEO for satellite system
*                           add BDS L1 code for RINEX 3.02 and RTCM 3.2
*                           support BDS L1 in satwavelen()
*           2014/05/29 1.27 fix bug on obs2code() to search obs code table
*           2014/08/26 1.28 fix problem on output of uncompress() for tar file
*                           add function to swap trace file with keywords
*           2014/10/21 1.29 strtok() -> strtok_r() in expath() for thread-safe
*                           add bdsmodear in procopt_default
*           2015/03/19 1.30 fix bug on interpolation of erp values in geterp()
*                           add leap second insertion before 2015/07/01 00:00
*                           add api read_leaps()
*           2015/05/31 1.31 delete api windupcorr()
*           2015/08/08 1.32 add compile option CPUTIME_IN_GPST
*                           add api add_fatal()
*                           support usno leapsec.dat for api read_leaps()
*           2016/01/23 1.33 enable septentrio
*           2016/02/05 1.34 support GLONASS for savenav(), loadnav()
*           2016/06/11 1.35 delete trace() in reppath() to avoid deadlock
*           2016/07/01 1.36 support IRNSS
*                           add leap second before 2017/1/1 00:00:00
*           2016/07/29 1.37 rename api compress() -> rtk_uncompress()
*                           rename api crc16()    -> rtk_crc16()
*                           rename api crc24q()   -> rtk_crc24q()
*                           rename api crc32()    -> rtk_crc32()
*           2016/08/20 1.38 fix type incompatibility in win64 environment
*                           change constant _POSIX_C_SOURCE 199309 -> 199506
*           2016/08/21 1.39 fix bug on week overflow in time2gpst()/gpst2time()
*           2016/09/05 1.40 fix bug on invalid nav data read in readnav()
*           2016/09/17 1.41 suppress warnings
*           2016/09/19 1.42 modify api deg2dms() to consider numerical error
*           2017/04/11 1.43 delete EXPORT for global variables
*           2018/10/10 1.44 modify api satexclude()
*           2020/11/30 1.45 add API code2idx() to get freq-index
*                           add API code2freq() to get carrier frequency
*                           add API timereset() to reset current time
*                           modify API obs2code(), code2obs() and setcodepri()
*                           delete API satwavelen()
*                           delete API csmooth()
*                           delete global variable lam_carr[]
*                           compensate L3,L4,... PCVs by L2 PCV if no PCV data
*                            in input file by API readpcv()
*                           add support hatanaka-compressed RINEX files with
*                            extension ".crx" or ".CRX"
*                           update stream format strings table
*                           update obs code strings and priority table
*                           use integer types in stdint.h
*                           surppress warnings
*-----------------------------------------------------------------------------*/
#define _POSIX_C_SOURCE 199506
#include <stdarg.h>
#include <ctype.h>
#include <errno.h>
#ifndef WIN32
#include <dirent.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

/* constants -----------------------------------------------------------------*/

#define POLYCRC32   0xEDB88320u /* CRC32 polynomial */
#define POLYCRC24Q  0x1864CFBu  /* CRC24Q polynomial */

#define SQR(x)      ((x)*(x))
#define MAX_VAR_EPH SQR(300.0)  /* max variance eph to reject satellite (m^2) */

static const double od_rtk_rtkcmn_gpst0[]={1980,1, 6,0,0,0}; /* gps time reference */
static const double od_rtk_rtkcmn_gst0 []={1999,8,22,0,0,0}; /* galileo system time reference */
static const double od_rtk_rtkcmn_bdt0 []={2006,1, 1,0,0,0}; /* beidou time reference */

static double od_rtk_rtkcmn_leaps[MAXLEAPS+1][7]={ /* leap seconds (y,m,d,h,m,s,utc-gpst) */
    {2017,1,1,0,0,0,-18},
    {2015,7,1,0,0,0,-17},
    {2012,7,1,0,0,0,-16},
    {2009,1,1,0,0,0,-15},
    {2006,1,1,0,0,0,-14},
    {1999,1,1,0,0,0,-13},
    {1997,7,1,0,0,0,-12},
    {1996,1,1,0,0,0,-11},
    {1994,7,1,0,0,0,-10},
    {1993,7,1,0,0,0, -9},
    {1992,7,1,0,0,0, -8},
    {1991,1,1,0,0,0, -7},
    {1990,1,1,0,0,0, -6},
    {1988,1,1,0,0,0, -5},
    {1985,7,1,0,0,0, -4},
    {1983,7,1,0,0,0, -3},
    {1982,7,1,0,0,0, -2},
    {1981,7,1,0,0,0, -1},
    {0}
};
const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
    10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
    31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
    46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
    61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
    74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
    88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
    101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
    113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
    126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
    138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
};
const prcopt_t prcopt_default={ /* defaults processing options */
    PMODE_SINGLE,0,2,SYS_GPS,   /* mode,soltype,nf,navsys */
    15.0*D2R,{{0,0}},           /* elmin,snrmask */
    0,1,1,1,                    /* sateph,modear,glomodear,bdsmodear */
    5,0,10,1,                   /* maxout,minlock,minfix,armaxiter */
    0,0,0,0,                    /* estion,esttrop,dynamics,tidecorr */
    1,0,0,0,0,                  /* niter,codesmooth,intpref,sbascorr,sbassatsel */
    0,0,                        /* rovpos,refpos */
    {100.0,100.0},              /* eratio[] */
    {100.0,0.003,0.003,0.0,1.0}, /* err[] */
    {30.0,0.03,0.3},            /* std[] */
    {1E-4,1E-3,1E-4,1E-1,1E-2,0.0}, /* prn[] */
    5E-12,                      /* sclkstab */
    {3.0,0.9999,0.25,0.1,0.05}, /* thresar */
    0.0,0.0,0.05,               /* elmaskar,almaskhold,thresslip */
    30.0,30.0,30.0,             /* maxtdif,maxinno,maxgdop */
    {0},{0},{0},                /* baseline,ru,rb */
    {"",""},                    /* anttype */
    {{0}},{{0}},{0}             /* antdel,pcv,exsats */
};
const solopt_t solopt_default={ /* defaults solution output options */
    SOLF_LLH,TIMES_GPST,1,3,    /* posf,times,timef,timeu */
    0,1,0,0,0,0,0,              /* degf,outhead,outopt,outvel,datum,height,geoid */
    0,0,0,                      /* solstatic,sstat,trace */
    {0.0,0.0},                  /* nmeaintv */
    " ",""                      /* separator/program name */
};
const char *formatstrs[32]={    /* stream format strings */
    "RTCM 2",                   /*  0 */
    "RTCM 3",                   /*  1 */
    "NovAtel OEM7",             /*  2 */
    "NovAtel OEM3",             /*  3 */
    "u-blox UBX",               /*  4 */
    "Superstar II",             /*  5 */
    "Hemisphere",               /*  6 */
    "SkyTraq",                  /*  7 */
    "Javad GREIS",              /*  8 */
    "NVS BINR",                 /*  9 */
    "BINEX",                    /* 10 */
    "Trimble RT17",             /* 11 */
    "Septentrio SBF",           /* 12 */
    "RINEX",                    /* 13 */
    "SP3",                      /* 14 */
    "RINEX CLK",                /* 15 */
    "SBAS",                     /* 16 */
    "NMEA 0183",                /* 17 */
    NULL
};
static char *od_rtk_rtkcmn_obscodes[]={       /* observation code strings */
    
    ""  ,"1C","1P","1W","1Y", "1M","1N","1S","1L","1E", /*  0- 9 */
    "1A","1B","1X","1Z","2C", "2D","2S","2L","2X","2P", /* 10-19 */
    "2W","2Y","2M","2N","5I", "5Q","5X","7I","7Q","7X", /* 20-29 */
    "6A","6B","6C","6X","6Z", "6S","6L","8L","8Q","8X", /* 30-39 */
    "2I","2Q","6I","6Q","3I", "3Q","3X","1I","1Q","5A", /* 40-49 */
    "5B","5C","9A","9B","9C", "9X","1D","5D","5P","5Z", /* 50-59 */
    "6E","7D","7P","7Z","8D", "8P","4A","4B","4X",""    /* 60-69 */
};
static char od_rtk_rtkcmn_codepris[7][MAXFREQ][16]={  /* code priority for each freq-index */
   /*    0         1          2          3         4         5     */
    {"CPYWMNSL","PYWCMNDLSX","IQX"     ,""       ,""       ,""      ,""}, /* GPS */
    {"CPABX"   ,"PCABX"     ,"IQX"     ,""       ,""       ,""      ,""}, /* GLO */
    {"CABXZ"   ,"IQX"       ,"IQX"     ,"ABCXZ"  ,"IQX"    ,""      ,""}, /* GAL */
    {"CLSXZ"   ,"LSX"       ,"IQXDPZ"  ,"LSXEZ"  ,""       ,""      ,""}, /* QZS */
    {"C"       ,"IQX"       ,""        ,""       ,""       ,""      ,""}, /* SBS */
    {"IQXDPAN" ,"IQXDPZ"    ,"DPX"     ,"IQXA"   ,"DPX"    ,""      ,""}, /* BDS */
    {"ABCX"    ,"ABCX"      ,""        ,""       ,""       ,""      ,""}  /* IRN */
};
static fatalfunc_t *od_rtk_rtkcmn_fatalfunc=NULL; /* fatal callback function */

/* crc tables generated by util/gencrc ---------------------------------------*/
static const uint16_t od_rtk_rtkcmn_tbl_CRC16[]={
    0x0000,0x1021,0x2042,0x3063,0x4084,0x50A5,0x60C6,0x70E7,
    0x8108,0x9129,0xA14A,0xB16B,0xC18C,0xD1AD,0xE1CE,0xF1EF,
    0x1231,0x0210,0x3273,0x2252,0x52B5,0x4294,0x72F7,0x62D6,
    0x9339,0x8318,0xB37B,0xA35A,0xD3BD,0xC39C,0xF3FF,0xE3DE,
    0x2462,0x3443,0x0420,0x1401,0x64E6,0x74C7,0x44A4,0x5485,
    0xA56A,0xB54B,0x8528,0x9509,0xE5EE,0xF5CF,0xC5AC,0xD58D,
    0x3653,0x2672,0x1611,0x0630,0x76D7,0x66F6,0x5695,0x46B4,
    0xB75B,0xA77A,0x9719,0x8738,0xF7DF,0xE7FE,0xD79D,0xC7BC,
    0x48C4,0x58E5,0x6886,0x78A7,0x0840,0x1861,0x2802,0x3823,
    0xC9CC,0xD9ED,0xE98E,0xF9AF,0x8948,0x9969,0xA90A,0xB92B,
    0x5AF5,0x4AD4,0x7AB7,0x6A96,0x1A71,0x0A50,0x3A33,0x2A12,
    0xDBFD,0xCBDC,0xFBBF,0xEB9E,0x9B79,0x8B58,0xBB3B,0xAB1A,
    0x6CA6,0x7C87,0x4CE4,0x5CC5,0x2C22,0x3C03,0x0C60,0x1C41,
    0xEDAE,0xFD8F,0xCDEC,0xDDCD,0xAD2A,0xBD0B,0x8D68,0x9D49,
    0x7E97,0x6EB6,0x5ED5,0x4EF4,0x3E13,0x2E32,0x1E51,0x0E70,
    0xFF9F,0xEFBE,0xDFDD,0xCFFC,0xBF1B,0xAF3A,0x9F59,0x8F78,
    0x9188,0x81A9,0xB1CA,0xA1EB,0xD10C,0xC12D,0xF14E,0xE16F,
    0x1080,0x00A1,0x30C2,0x20E3,0x5004,0x4025,0x7046,0x6067,
    0x83B9,0x9398,0xA3FB,0xB3DA,0xC33D,0xD31C,0xE37F,0xF35E,
    0x02B1,0x1290,0x22F3,0x32D2,0x4235,0x5214,0x6277,0x7256,
    0xB5EA,0xA5CB,0x95A8,0x8589,0xF56E,0xE54F,0xD52C,0xC50D,
    0x34E2,0x24C3,0x14A0,0x0481,0x7466,0x6447,0x5424,0x4405,
    0xA7DB,0xB7FA,0x8799,0x97B8,0xE75F,0xF77E,0xC71D,0xD73C,
    0x26D3,0x36F2,0x0691,0x16B0,0x6657,0x7676,0x4615,0x5634,
    0xD94C,0xC96D,0xF90E,0xE92F,0x99C8,0x89E9,0xB98A,0xA9AB,
    0x5844,0x4865,0x7806,0x6827,0x18C0,0x08E1,0x3882,0x28A3,
    0xCB7D,0xDB5C,0xEB3F,0xFB1E,0x8BF9,0x9BD8,0xABBB,0xBB9A,
    0x4A75,0x5A54,0x6A37,0x7A16,0x0AF1,0x1AD0,0x2AB3,0x3A92,
    0xFD2E,0xED0F,0xDD6C,0xCD4D,0xBDAA,0xAD8B,0x9DE8,0x8DC9,
    0x7C26,0x6C07,0x5C64,0x4C45,0x3CA2,0x2C83,0x1CE0,0x0CC1,
    0xEF1F,0xFF3E,0xCF5D,0xDF7C,0xAF9B,0xBFBA,0x8FD9,0x9FF8,
    0x6E17,0x7E36,0x4E55,0x5E74,0x2E93,0x3EB2,0x0ED1,0x1EF0
};
static const uint32_t od_rtk_rtkcmn_tbl_CRC24Q[]={
    0x000000,0x864CFB,0x8AD50D,0x0C99F6,0x93E6E1,0x15AA1A,0x1933EC,0x9F7F17,
    0xA18139,0x27CDC2,0x2B5434,0xAD18CF,0x3267D8,0xB42B23,0xB8B2D5,0x3EFE2E,
    0xC54E89,0x430272,0x4F9B84,0xC9D77F,0x56A868,0xD0E493,0xDC7D65,0x5A319E,
    0x64CFB0,0xE2834B,0xEE1ABD,0x685646,0xF72951,0x7165AA,0x7DFC5C,0xFBB0A7,
    0x0CD1E9,0x8A9D12,0x8604E4,0x00481F,0x9F3708,0x197BF3,0x15E205,0x93AEFE,
    0xAD50D0,0x2B1C2B,0x2785DD,0xA1C926,0x3EB631,0xB8FACA,0xB4633C,0x322FC7,
    0xC99F60,0x4FD39B,0x434A6D,0xC50696,0x5A7981,0xDC357A,0xD0AC8C,0x56E077,
    0x681E59,0xEE52A2,0xE2CB54,0x6487AF,0xFBF8B8,0x7DB443,0x712DB5,0xF7614E,
    0x19A3D2,0x9FEF29,0x9376DF,0x153A24,0x8A4533,0x0C09C8,0x00903E,0x86DCC5,
    0xB822EB,0x3E6E10,0x32F7E6,0xB4BB1D,0x2BC40A,0xAD88F1,0xA11107,0x275DFC,
    0xDCED5B,0x5AA1A0,0x563856,0xD074AD,0x4F0BBA,0xC94741,0xC5DEB7,0x43924C,
    0x7D6C62,0xFB2099,0xF7B96F,0x71F594,0xEE8A83,0x68C678,0x645F8E,0xE21375,
    0x15723B,0x933EC0,0x9FA736,0x19EBCD,0x8694DA,0x00D821,0x0C41D7,0x8A0D2C,
    0xB4F302,0x32BFF9,0x3E260F,0xB86AF4,0x2715E3,0xA15918,0xADC0EE,0x2B8C15,
    0xD03CB2,0x567049,0x5AE9BF,0xDCA544,0x43DA53,0xC596A8,0xC90F5E,0x4F43A5,
    0x71BD8B,0xF7F170,0xFB6886,0x7D247D,0xE25B6A,0x641791,0x688E67,0xEEC29C,
    0x3347A4,0xB50B5F,0xB992A9,0x3FDE52,0xA0A145,0x26EDBE,0x2A7448,0xAC38B3,
    0x92C69D,0x148A66,0x181390,0x9E5F6B,0x01207C,0x876C87,0x8BF571,0x0DB98A,
    0xF6092D,0x7045D6,0x7CDC20,0xFA90DB,0x65EFCC,0xE3A337,0xEF3AC1,0x69763A,
    0x578814,0xD1C4EF,0xDD5D19,0x5B11E2,0xC46EF5,0x42220E,0x4EBBF8,0xC8F703,
    0x3F964D,0xB9DAB6,0xB54340,0x330FBB,0xAC70AC,0x2A3C57,0x26A5A1,0xA0E95A,
    0x9E1774,0x185B8F,0x14C279,0x928E82,0x0DF195,0x8BBD6E,0x872498,0x016863,
    0xFAD8C4,0x7C943F,0x700DC9,0xF64132,0x693E25,0xEF72DE,0xE3EB28,0x65A7D3,
    0x5B59FD,0xDD1506,0xD18CF0,0x57C00B,0xC8BF1C,0x4EF3E7,0x426A11,0xC426EA,
    0x2AE476,0xACA88D,0xA0317B,0x267D80,0xB90297,0x3F4E6C,0x33D79A,0xB59B61,
    0x8B654F,0x0D29B4,0x01B042,0x87FCB9,0x1883AE,0x9ECF55,0x9256A3,0x141A58,
    0xEFAAFF,0x69E604,0x657FF2,0xE33309,0x7C4C1E,0xFA00E5,0xF69913,0x70D5E8,
    0x4E2BC6,0xC8673D,0xC4FECB,0x42B230,0xDDCD27,0x5B81DC,0x57182A,0xD154D1,
    0x26359F,0xA07964,0xACE092,0x2AAC69,0xB5D37E,0x339F85,0x3F0673,0xB94A88,
    0x87B4A6,0x01F85D,0x0D61AB,0x8B2D50,0x145247,0x921EBC,0x9E874A,0x18CBB1,
    0xE37B16,0x6537ED,0x69AE1B,0xEFE2E0,0x709DF7,0xF6D10C,0xFA48FA,0x7C0401,
    0x42FA2F,0xC4B6D4,0xC82F22,0x4E63D9,0xD11CCE,0x575035,0x5BC9C3,0xDD8538
};
/* function prototypes -------------------------------------------------------*/
#ifdef MKL
#define LAPACK
#define dgemm_      dgemm
#define dgetrf_     dgetrf
#define dgetri_     dgetri
#define dgetrs_     dgetrs
#endif
#ifdef LAPACK
extern void dgemm_(char *, char *, int *, int *, int *, double *, double *,
                   int *, double *, int *, double *, double *, int *);
extern void dgetrf_(int *, int *, double *, int *, int *, int *);
extern void dgetri_(int *, double *, int *, int *, double *, int *, int *);
extern void dgetrs_(char *, int *, int *, double *, int *, int *, double *,
                    int *, int *);
#endif

#ifdef IERS_MODEL
extern int gmf_(double *mjd, double *lat, double *lon, double *hgt, double *zd,
                double *gmfh, double *gmfw);
#endif

/* fatal error ---------------------------------------------------------------*/
static void od_rtk_rtkcmn_fatalerr(const char *format, ...)
{
    char msg[1024];
    va_list ap;
    va_start(ap,format); vsprintf(msg,format,ap); va_end(ap);
    if (od_rtk_rtkcmn_fatalfunc) od_rtk_rtkcmn_fatalfunc(msg);
    else fprintf(stderr,"%s",msg);
    exit(-9);
}
/* add fatal callback function -------------------------------------------------
* add fatal callback function for mat(),zeros(),imat()
* args   : fatalfunc_t *func I  callback function
* return : none
* notes  : if malloc() failed in return : none
*-----------------------------------------------------------------------------*/
extern void add_fatal(fatalfunc_t *func)
{
    od_rtk_rtkcmn_fatalfunc=func;
}
/* satellite system+prn/slot number to satellite number ------------------------
* convert satellite system+prn/slot number to satellite number
* args   : int    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
*          int    prn       I   satellite prn/slot number
* return : satellite number (0:error)
*-----------------------------------------------------------------------------*/
extern int satno(int sys, int prn)
{
    if (prn<=0) return 0;
    switch (sys) {
        case SYS_GPS:
            if (prn<MINPRNGPS||MAXPRNGPS<prn) return 0;
            return prn-MINPRNGPS+1;
        case SYS_GLO:
            if (prn<MINPRNGLO||MAXPRNGLO<prn) return 0;
            return NSATGPS+prn-MINPRNGLO+1;
        case SYS_GAL:
            if (prn<MINPRNGAL||MAXPRNGAL<prn) return 0;
            return NSATGPS+NSATGLO+prn-MINPRNGAL+1;
        case SYS_QZS:
            if (prn<MINPRNQZS||MAXPRNQZS<prn) return 0;
            return NSATGPS+NSATGLO+NSATGAL+prn-MINPRNQZS+1;
        case SYS_CMP:
            if (prn<MINPRNCMP||MAXPRNCMP<prn) return 0;
            return NSATGPS+NSATGLO+NSATGAL+NSATQZS+prn-MINPRNCMP+1;
        case SYS_IRN:
            if (prn<MINPRNIRN||MAXPRNIRN<prn) return 0;
            return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+prn-MINPRNIRN+1;
        case SYS_LEO:
            if (prn<MINPRNLEO||MAXPRNLEO<prn) return 0;
            return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+
                   prn-MINPRNLEO+1;
        case SYS_SBS:
            if (prn<MINPRNSBS||MAXPRNSBS<prn) return 0;
            return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATLEO+
                   prn-MINPRNSBS+1;
    }
    return 0;
}
/* satellite number to satellite system ----------------------------------------
* convert satellite number to satellite system
* args   : int    sat       I   satellite number (1-MAXSAT)
*          int    *prn      IO  satellite prn/slot number (NULL: no output)
* return : satellite system (SYS_GPS,SYS_GLO,...)
*-----------------------------------------------------------------------------*/
extern int satsys(int sat, int *prn)
{
    int sys=SYS_NONE;
    if (sat<=0||MAXSAT<sat) sat=0;
    else if (sat<=NSATGPS) {
        sys=SYS_GPS; sat+=MINPRNGPS-1;
    }
    else if ((sat-=NSATGPS)<=NSATGLO) {
        sys=SYS_GLO; sat+=MINPRNGLO-1;
    }
    else if ((sat-=NSATGLO)<=NSATGAL) {
        sys=SYS_GAL; sat+=MINPRNGAL-1;
    }
    else if ((sat-=NSATGAL)<=NSATQZS) {
        sys=SYS_QZS; sat+=MINPRNQZS-1; 
    }
    else if ((sat-=NSATQZS)<=NSATCMP) {
        sys=SYS_CMP; sat+=MINPRNCMP-1; 
    }
    else if ((sat-=NSATCMP)<=NSATIRN) {
        sys=SYS_IRN; sat+=MINPRNIRN-1; 
    }
    else if ((sat-=NSATIRN)<=NSATLEO) {
        sys=SYS_LEO; sat+=MINPRNLEO-1; 
    }
    else if ((sat-=NSATLEO)<=NSATSBS) {
        sys=SYS_SBS; sat+=MINPRNSBS-1; 
    }
    else sat=0;
    if (prn) *prn=sat;
    return sys;
}
/* satellite id to satellite number --------------------------------------------
* convert satellite id to satellite number
* args   : char   *id       I   satellite id (nn,Gnn,Rnn,Enn,Jnn,Cnn,Inn or Snn)
* return : satellite number (0: error)
* notes  : 120-142 and 193-199 are also recognized as sbas and qzss
*-----------------------------------------------------------------------------*/
extern int satid2no(const char *id)
{
    int sys,prn;
    char code;
    
    if (sscanf(id,"%d",&prn)==1) {
        if      (MINPRNGPS<=prn&&prn<=MAXPRNGPS) sys=SYS_GPS;
        else if (MINPRNSBS<=prn&&prn<=MAXPRNSBS) sys=SYS_SBS;
        else if (MINPRNQZS<=prn&&prn<=MAXPRNQZS) sys=SYS_QZS;
        else return 0;
        return satno(sys,prn);
    }
    if (sscanf(id,"%c%d",&code,&prn)<2) return 0;
    
    switch (code) {
        case 'G': sys=SYS_GPS; prn+=MINPRNGPS-1; break;
        case 'R': sys=SYS_GLO; prn+=MINPRNGLO-1; break;
        case 'E': sys=SYS_GAL; prn+=MINPRNGAL-1; break;
        case 'J': sys=SYS_QZS; prn+=MINPRNQZS-1; break;
        case 'C': sys=SYS_CMP; prn+=MINPRNCMP-1; break;
        case 'I': sys=SYS_IRN; prn+=MINPRNIRN-1; break;
        case 'L': sys=SYS_LEO; prn+=MINPRNLEO-1; break;
        case 'S': sys=SYS_SBS; prn+=100; break;
        default: return 0;
    }
    return satno(sys,prn);
}
/* satellite number to satellite id --------------------------------------------
* convert satellite number to satellite id
* args   : int    sat       I   satellite number
*          char   *id       O   satellite id (Gnn,Rnn,Enn,Jnn,Cnn,Inn or nnn)
* return : none
*-----------------------------------------------------------------------------*/
extern void satno2id(int sat, char *id)
{
    int prn;
    switch (satsys(sat,&prn)) {
        case SYS_GPS: sprintf(id,"G%02d",prn-MINPRNGPS+1); return;
        case SYS_GLO: sprintf(id,"R%02d",prn-MINPRNGLO+1); return;
        case SYS_GAL: sprintf(id,"E%02d",prn-MINPRNGAL+1); return;
        case SYS_QZS: sprintf(id,"J%02d",prn-MINPRNQZS+1); return;
        case SYS_CMP: sprintf(id,"C%02d",prn-MINPRNCMP+1); return;
        case SYS_IRN: sprintf(id,"I%02d",prn-MINPRNIRN+1); return;
        case SYS_LEO: sprintf(id,"L%02d",prn-MINPRNLEO+1); return;
        case SYS_SBS: sprintf(id,"%03d" ,prn); return;
    }
    strcpy(id,"");
}
/* test excluded satellite -----------------------------------------------------
* test excluded satellite
* args   : int    sat       I   satellite number
*          double var       I   variance of ephemeris (m^2)
*          int    svh       I   sv health flag
*          prcopt_t *opt    I   processing options (NULL: not used)
* return : status (1:excluded,0:not excluded)
*-----------------------------------------------------------------------------*/
extern int satexclude(int sat, double var, int svh, const prcopt_t *opt)
{
    int sys=satsys(sat,NULL);
    
    if (svh<0) return 1; /* ephemeris unavailable */
    
    if (opt) {
        if (opt->exsats[sat-1]==1) return 1; /* excluded satellite */
        if (opt->exsats[sat-1]==2) return 0; /* included satellite */
        if (!(sys&opt->navsys)) return 1; /* unselected sat sys */
    }
    if (sys==SYS_QZS) svh&=0xFE; /* mask QZSS LEX health */
    if (svh) {
        trace(3,"unhealthy satellite: sat=%3d svh=%02X\n",sat,svh);
        return 1;
    }
    if (var>MAX_VAR_EPH) {
        trace(3,"invalid ura satellite: sat=%3d ura=%.2f\n",sat,sqrt(var));
        return 1;
    }
    return 0;
}
/* test SNR mask ---------------------------------------------------------------
* test SNR mask
* args   : int    base      I   rover or base-station (0:rover,1:base station)
*          int    idx       I   frequency index (0:L1,1:L2,2:L3,...)
*          double el        I   elevation angle (rad)
*          double snr       I   C/N0 (dBHz)
*          snrmask_t *mask  I   SNR mask
* return : status (1:masked,0:unmasked)
*-----------------------------------------------------------------------------*/
extern int testsnr(int base, int idx, double el, double snr,
                   const snrmask_t *mask)
{
    double minsnr,a;
    int i;
    
    if (!mask->ena[base]||idx<0||idx>=NFREQ) return 0;
    
    a=(el*R2D+5.0)/10.0;
    i=(int)floor(a); a-=i;
    if      (i<1) minsnr=mask->mask[idx][0];
    else if (i>8) minsnr=mask->mask[idx][8];
    else minsnr=(1.0-a)*mask->mask[idx][i-1]+a*mask->mask[idx][i];
    
    return snr<minsnr;
}
/* obs type string to obs code -------------------------------------------------
* convert obs code type string to obs code
* args   : char   *str      I   obs code string ("1C","1P","1Y",...)
* return : obs code (CODE_???)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
extern uint8_t obs2code(const char *obs)
{
    int i;
    
    for (i=1;*od_rtk_rtkcmn_obscodes[i];i++) {
        if (strcmp(od_rtk_rtkcmn_obscodes[i],obs)) continue;
        return (uint8_t)i;
    }
    return CODE_NONE;
}
/* obs code to obs code string -------------------------------------------------
* convert obs code to obs code string
* args   : uint8_t code     I   obs code (CODE_???)
* return : obs code string ("1C","1P","1P",...)
* notes  : obs codes are based on RINEX 3.04
*-----------------------------------------------------------------------------*/
extern char *code2obs(uint8_t code)
{
    if (code<=CODE_NONE||MAXCODE<code) return "";
    return od_rtk_rtkcmn_obscodes[code];
}
/* GPS obs code to frequency -------------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_GPS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '2': *freq=FREQ2; return 1; /* L2 */
        case '5': *freq=FREQ5; return 2; /* L5 */
    }
    return -1;
}
/* GLONASS obs code to frequency ---------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_GLO(uint8_t code, int fcn, double *freq)
{
    char *obs=code2obs(code);
    
    if (fcn<-7||fcn>6) return -1;
    
    switch (obs[0]) {
        case '1': *freq=FREQ1_GLO+DFRQ1_GLO*fcn; return 0; /* G1 */
        case '2': *freq=FREQ2_GLO+DFRQ2_GLO*fcn; return 1; /* G2 */
        case '3': *freq=FREQ3_GLO;               return 2; /* G3 */
        case '4': *freq=FREQ1a_GLO;              return 0; /* G1a */
        case '6': *freq=FREQ2a_GLO;              return 1; /* G2a */
    }
    return -1;
}
/* Galileo obs code to frequency ---------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_GAL(uint8_t code, double *freq)
{
    char *obs=code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* E1 */
        case '7': *freq=FREQ7; return 1; /* E5b */
        case '5': *freq=FREQ5; return 2; /* E5a */
        case '6': *freq=FREQ6; return 3; /* E6 */
        case '8': *freq=FREQ8; return 4; /* E5ab */
    }
    return -1;
}
/* QZSS obs code to frequency ------------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_QZS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '2': *freq=FREQ2; return 1; /* L2 */
        case '5': *freq=FREQ5; return 2; /* L5 */
        case '6': *freq=FREQ6; return 3; /* L6 */
    }
    return -1;
}
/* SBAS obs code to frequency ------------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_SBS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1; return 0; /* L1 */
        case '5': *freq=FREQ5; return 1; /* L5 */
    }
    return -1;
}
/* BDS obs code to frequency -------------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_BDS(uint8_t code, double *freq)
{
    char *obs=code2obs(code);
    
    switch (obs[0]) {
        case '1': *freq=FREQ1;     return 0; /* B1C */
        case '2': *freq=FREQ1_CMP; return 0; /* B1I */
        case '7': *freq=FREQ2_CMP; return 1; /* B2I/B2b */
        case '5': *freq=FREQ5;     return 2; /* B2a */
        case '6': *freq=FREQ3_CMP; return 3; /* B3 */
        case '8': *freq=FREQ8;     return 4; /* B2ab */
    }
    return -1;
}
/* NavIC obs code to frequency -----------------------------------------------*/
static int od_rtk_rtkcmn_code2freq_IRN(uint8_t code, double *freq)
{
    char *obs=code2obs(code);
    
    switch (obs[0]) {
        case '5': *freq=FREQ5; return 0; /* L5 */
        case '9': *freq=FREQ9; return 1; /* S */
    }
    return -1;
}
/* system and obs code to frequency index --------------------------------------
* convert system and obs code to frequency index
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
* return : frequency index (-1: error)
*                       0     1     2     3     4 
*           --------------------------------------
*            GPS       L1    L2    L5     -     - 
*            GLONASS   G1    G2    G3     -     -  (G1=G1,G1a,G2=G2,G2a)
*            Galileo   E1    E5b   E5a   E6   E5ab
*            QZSS      L1    L2    L5    L6     - 
*            SBAS      L1     -    L5     -     -
*            BDS       B1    B2    B2a   B3   B2ab (B1=B1I,B1C,B2=B2I,B2b)
*            NavIC     L5     S     -     -     - 
*-----------------------------------------------------------------------------*/
extern int code2idx(int sys, uint8_t code)
{
    double freq;
    
    switch (sys) {
        case SYS_GPS: return od_rtk_rtkcmn_code2freq_GPS(code,&freq);
        case SYS_GLO: return od_rtk_rtkcmn_code2freq_GLO(code,0,&freq);
        case SYS_GAL: return od_rtk_rtkcmn_code2freq_GAL(code,&freq);
        case SYS_QZS: return od_rtk_rtkcmn_code2freq_QZS(code,&freq);
        case SYS_SBS: return od_rtk_rtkcmn_code2freq_SBS(code,&freq);
        case SYS_CMP: return od_rtk_rtkcmn_code2freq_BDS(code,&freq);
        case SYS_IRN: return od_rtk_rtkcmn_code2freq_IRN(code,&freq);
    }
    return -1;
}
/* system and obs code to frequency --------------------------------------------
* convert system and obs code to carrier frequency
* args   : int    sys       I   satellite system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
*          int    fcn       I   frequency channel number for GLONASS
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
extern double code2freq(int sys, uint8_t code, int fcn)
{
    double freq=0.0;
    
    switch (sys) {
        case SYS_GPS: (void)od_rtk_rtkcmn_code2freq_GPS(code,&freq); break;
        case SYS_GLO: (void)od_rtk_rtkcmn_code2freq_GLO(code,fcn,&freq); break;
        case SYS_GAL: (void)od_rtk_rtkcmn_code2freq_GAL(code,&freq); break;
        case SYS_QZS: (void)od_rtk_rtkcmn_code2freq_QZS(code,&freq); break;
        case SYS_SBS: (void)od_rtk_rtkcmn_code2freq_SBS(code,&freq); break;
        case SYS_CMP: (void)od_rtk_rtkcmn_code2freq_BDS(code,&freq); break;
        case SYS_IRN: (void)od_rtk_rtkcmn_code2freq_IRN(code,&freq); break;
    }
    return freq;
}
/* satellite and obs code to frequency -----------------------------------------
* convert satellite and obs code to carrier frequency
* args   : int    sat       I   satellite number
*          uint8_t code     I   obs code (CODE_???)
*          nav_t  *nav_t    I   navigation data for GLONASS (NULL: not used)
* return : carrier frequency (Hz) (0.0: error)
*-----------------------------------------------------------------------------*/
extern double sat2freq(int sat, uint8_t code, const nav_t *nav)
{
    int i,fcn=0,sys,prn;
    
    sys=satsys(sat,&prn);
    
    if (sys==SYS_GLO) {
        if (!nav) return 0.0;
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        if (i<nav->ng) {
            fcn=nav->geph[i].frq;
        }
        else if (nav->glo_fcn[prn-1]>0) {
            fcn=nav->glo_fcn[prn-1]-8;
        }
        else return 0.0;
    }
    return code2freq(sys,code,fcn);
}
/* set code priority -----------------------------------------------------------
* set code priority for multiple codes in a frequency
* args   : int    sys       I   system (or of SYS_???)
*          int    idx       I   frequency index (0- )
*          char   *pri      I   priority of codes (series of code characters)
*                               (higher priority precedes lower)
* return : none
*-----------------------------------------------------------------------------*/
extern void setcodepri(int sys, int idx, const char *pri)
{
    trace(3,"setcodepri:sys=%d idx=%d pri=%s\n",sys,idx,pri);
    
    if (idx<0||idx>=MAXFREQ) return;
    if (sys&SYS_GPS) strcpy(od_rtk_rtkcmn_codepris[0][idx],pri);
    if (sys&SYS_GLO) strcpy(od_rtk_rtkcmn_codepris[1][idx],pri);
    if (sys&SYS_GAL) strcpy(od_rtk_rtkcmn_codepris[2][idx],pri);
    if (sys&SYS_QZS) strcpy(od_rtk_rtkcmn_codepris[3][idx],pri);
    if (sys&SYS_SBS) strcpy(od_rtk_rtkcmn_codepris[4][idx],pri);
    if (sys&SYS_CMP) strcpy(od_rtk_rtkcmn_codepris[5][idx],pri);
    if (sys&SYS_IRN) strcpy(od_rtk_rtkcmn_codepris[6][idx],pri);
}
/* get code priority -----------------------------------------------------------
* get code priority for multiple codes in a frequency
* args   : int    sys       I   system (SYS_???)
*          uint8_t code     I   obs code (CODE_???)
*          char   *opt      I   code options (NULL:no option)
* return : priority (15:highest-1:lowest,0:error)
*-----------------------------------------------------------------------------*/
extern int getcodepri(int sys, uint8_t code, const char *opt)
{
    const char *p,*optstr;
    char *obs,str[8]="";
    int i,j;
    
    switch (sys) {
        case SYS_GPS: i=0; optstr="-GL%2s"; break;
        case SYS_GLO: i=1; optstr="-RL%2s"; break;
        case SYS_GAL: i=2; optstr="-EL%2s"; break;
        case SYS_QZS: i=3; optstr="-JL%2s"; break;
        case SYS_SBS: i=4; optstr="-SL%2s"; break;
        case SYS_CMP: i=5; optstr="-CL%2s"; break;
        case SYS_IRN: i=6; optstr="-IL%2s"; break;
        default: return 0;
    }
    if ((j=code2idx(sys,code))<0) return 0;
    obs=code2obs(code);
    
    /* parse code options */
    for (p=opt;p&&(p=strchr(p,'-'));p++) {
        if (sscanf(p,optstr,str)<1||str[0]!=obs[0]) continue;
        return str[1]==obs[1]?15:0;
    }
    /* search code priority */
    return (p=strchr(od_rtk_rtkcmn_codepris[i][j],obs[1]))?14-(int)(p-od_rtk_rtkcmn_codepris[i][j]):0;
}
/* extract unsigned/signed bits ------------------------------------------------
* extract unsigned/signed bits from byte data
* args   : uint8_t *buff    I   byte data
*          int    pos       I   bit position from start of data (bits)
*          int    len       I   bit length (bits) (len<=32)
* return : extracted unsigned/signed bits
*-----------------------------------------------------------------------------*/
extern uint32_t getbitu(const uint8_t *buff, int pos, int len)
{
    uint32_t bits=0;
    int i;
    for (i=pos;i<pos+len;i++) bits=(bits<<1)+((buff[i/8]>>(7-i%8))&1u);
    return bits;
}
extern int32_t getbits(const uint8_t *buff, int pos, int len)
{
    uint32_t bits=getbitu(buff,pos,len);
    if (len<=0||32<=len||!(bits&(1u<<(len-1)))) return (int32_t)bits;
    return (int32_t)(bits|(~0u<<len)); /* extend sign */
}
/* set unsigned/signed bits ----------------------------------------------------
* set unsigned/signed bits to byte data
* args   : uint8_t *buff IO byte data
*          int    pos       I   bit position from start of data (bits)
*          int    len       I   bit length (bits) (len<=32)
*          [u]int32_t data  I   unsigned/signed data
* return : none
*-----------------------------------------------------------------------------*/
extern void setbitu(uint8_t *buff, int pos, int len, uint32_t data)
{
    uint32_t mask=1u<<(len-1);
    int i;
    if (len<=0||32<len) return;
    for (i=pos;i<pos+len;i++,mask>>=1) {
        if (data&mask) buff[i/8]|=1u<<(7-i%8); else buff[i/8]&=~(1u<<(7-i%8));
    }
}
extern void setbits(uint8_t *buff, int pos, int len, int32_t data)
{
    if (data<0) data|=1<<(len-1); else data&=~(1<<(len-1)); /* set sign bit */
    setbitu(buff,pos,len,(uint32_t)data);
}
/* crc-32 parity ---------------------------------------------------------------
* compute crc-32 parity for novatel raw
* args   : uint8_t *buff    I   data
*          int    len       I   data length (bytes)
* return : crc-32 parity
* notes  : see NovAtel OEMV firmware manual 1.7 32-bit CRC
*-----------------------------------------------------------------------------*/
extern uint32_t rtk_crc32(const uint8_t *buff, int len)
{
    uint32_t crc=0;
    int i,j;
    
    trace(4,"rtk_crc32: len=%d\n",len);
    
    for (i=0;i<len;i++) {
        crc^=buff[i];
        for (j=0;j<8;j++) {
            if (crc&1) crc=(crc>>1)^POLYCRC32; else crc>>=1;
        }
    }
    return crc;
}
/* crc-24q parity --------------------------------------------------------------
* compute crc-24q parity for sbas, rtcm3
* args   : uint8_t *buff    I   data
*          int    len       I   data length (bytes)
* return : crc-24Q parity
* notes  : see reference [2] A.4.3.3 Parity
*-----------------------------------------------------------------------------*/
extern uint32_t rtk_crc24q(const uint8_t *buff, int len)
{
    uint32_t crc=0;
    int i;
    
    trace(4,"rtk_crc24q: len=%d\n",len);
    
    for (i=0;i<len;i++) crc=((crc<<8)&0xFFFFFF)^od_rtk_rtkcmn_tbl_CRC24Q[(crc>>16)^buff[i]];
    return crc;
}
/* crc-16 parity ---------------------------------------------------------------
* compute crc-16 parity for binex, nvs
* args   : uint8_t *buff    I   data
*          int    len       I   data length (bytes)
* return : crc-16 parity
* notes  : see reference [10] A.3.
*-----------------------------------------------------------------------------*/
extern uint16_t rtk_crc16(const uint8_t *buff, int len)
{
    uint16_t crc=0;
    int i;
    
    trace(4,"rtk_crc16: len=%d\n",len);
    
    for (i=0;i<len;i++) {
        crc=(crc<<8)^od_rtk_rtkcmn_tbl_CRC16[((crc>>8)^buff[i])&0xFF];
    }
    return crc;
}
/* decode navigation data word -------------------------------------------------
* check party and decode navigation data word
* args   : uint32_t word    I   navigation data word (2+30bit)
*                               (previous word D29*-30* + current word D1-30)
*          uint8_t *data    O   decoded navigation data without parity
*                               (8bitx3)
* return : status (1:ok,0:parity error)
* notes  : see reference [1] 20.3.5.2 user parity algorithm
*-----------------------------------------------------------------------------*/
extern int decode_word(uint32_t word, uint8_t *data)
{
    const uint32_t hamming[]={
        0xBB1F3480,0x5D8F9A40,0xAEC7CD00,0x5763E680,0x6BB1F340,0x8B7A89C0
    };
    uint32_t parity=0,w;
    int i;
    
    trace(5,"decodeword: word=%08x\n",word);
    
    if (word&0x40000000) word^=0x3FFFFFC0;
    
    for (i=0;i<6;i++) {
        parity<<=1;
        for (w=(word&hamming[i])>>6;w;w>>=1) parity^=w&1;
    }
    if (parity!=(word&0x3F)) return 0;
    
    for (i=0;i<3;i++) data[i]=(uint8_t)(word>>(22-i*8));
    return 1;
}
/* new matrix ------------------------------------------------------------------
* allocate memory of matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *mat(int n, int m)
{
    double *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)malloc(sizeof(double)*n*m))) {
        od_rtk_rtkcmn_fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* new integer matrix ----------------------------------------------------------
* allocate memory of integer matrix 
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern int *imat(int n, int m)
{
    int *p;
    
    if (n<=0||m<=0) return NULL;
    if (!(p=(int *)malloc(sizeof(int)*n*m))) {
        od_rtk_rtkcmn_fatalerr("integer matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
    return p;
}
/* zero matrix -----------------------------------------------------------------
* generate new zero matrix
* args   : int    n,m       I   number of rows and columns of matrix
* return : matrix pointer (if n<=0 or m<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *zeros(int n, int m)
{
    double *p;
    
#if NOCALLOC
    if ((p=mat(n,m))) for (n=n*m-1;n>=0;n--) p[n]=0.0;
#else
    if (n<=0||m<=0) return NULL;
    if (!(p=(double *)calloc(n*m,sizeof(double)))) {
        od_rtk_rtkcmn_fatalerr("matrix memory allocation error: n=%d,m=%d\n",n,m);
    }
#endif
    return p;
}
/* identity matrix -------------------------------------------------------------
* generate new identity matrix
* args   : int    n         I   number of rows and columns of matrix
* return : matrix pointer (if n<=0, return NULL)
*-----------------------------------------------------------------------------*/
extern double *eye(int n)
{
    double *p;
    int i;
    
    if ((p=zeros(n,n))) for (i=0;i<n;i++) p[i+i*n]=1.0;
    return p;
}
/* inner product ---------------------------------------------------------------
* inner product of vectors
* args   : double *a,*b     I   vector a,b (n x 1)
*          int    n         I   size of vector a,b
* return : a'*b
*-----------------------------------------------------------------------------*/
extern double dot(const double *a, const double *b, int n)
{
    double c=0.0;
    
    while (--n>=0) c+=a[n]*b[n];
    return c;
}
/* euclid norm -----------------------------------------------------------------
* euclid norm of vector
* args   : double *a        I   vector a (n x 1)
*          int    n         I   size of vector a
* return : || a ||
*-----------------------------------------------------------------------------*/
extern double norm(const double *a, int n)
{
    return sqrt(dot(a,a,n));
}
/* outer product of 3d vectors -------------------------------------------------
* outer product of 3d vectors 
* args   : double *a,*b     I   vector a,b (3 x 1)
*          double *c        O   outer product (a x b) (3 x 1)
* return : none
*-----------------------------------------------------------------------------*/
extern void cross3(const double *a, const double *b, double *c)
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}
/* normalize 3d vector ---------------------------------------------------------
* normalize 3d vector
* args   : double *a        I   vector a (3 x 1)
*          double *b        O   normlized vector (3 x 1) || b || = 1
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int normv3(const double *a, double *b)
{
    double r;
    if ((r=norm(a,3))<=0.0) return 0;
    b[0]=a[0]/r;
    b[1]=a[1]/r;
    b[2]=a[2]/r;
    return 1;
}
/* copy matrix -----------------------------------------------------------------
* copy matrix
* args   : double *A        O   destination matrix A (n x m)
*          double *B        I   source matrix B (n x m)
*          int    n,m       I   number of rows and columns of matrix
* return : none
*-----------------------------------------------------------------------------*/
extern void matcpy(double *A, const double *B, int n, int m)
{
    memcpy(A,B,sizeof(double)*n*m);
}
/* matrix routines -----------------------------------------------------------*/

#ifdef LAPACK /* with LAPACK/BLAS or MKL */

/* multiply matrix (wrapper of blas dgemm) -------------------------------------
* multiply matrix by matrix (C=alpha*A*B+beta*C)
* args   : char   *tr       I  transpose flags ("N":normal,"T":transpose)
*          int    n,k,m     I  size of (transposed) matrix A,B
*          double alpha     I  alpha
*          double *A,*B     I  (transposed) matrix A (n x m), B (m x k)
*          double beta      I  beta
*          double *C        IO matrix C (n x k)
* return : none
*-----------------------------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    int lda=tr[0]=='T'?m:n,ldb=tr[1]=='T'?k:m;
    
    dgemm_((char *)tr,(char *)tr+1,&n,&k,&m,&alpha,(double *)A,&lda,(double *)B,
           &ldb,&beta,C,&n);
}
/* inverse of matrix -----------------------------------------------------------
* inverse of matrix (A=A^-1)
* args   : double *A        IO  matrix (n x n)
*          int    n         I   size of matrix A
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double *work;
    int info,lwork=n*16,*ipiv=imat(n,1);
    
    work=mat(lwork,1);
    dgetrf_(&n,&n,A,&n,ipiv,&info);
    if (!info) dgetri_(&n,A,&n,ipiv,work,&lwork,&info);
    free(ipiv); free(work);
    return info;
}
/* solve linear equation -------------------------------------------------------
* solve linear equation (X=A\Y or X=A'\Y)
* args   : char   *tr       I   transpose flag ("N":normal,"T":transpose)
*          double *A        I   input matrix A (n x n)
*          double *Y        I   input matrix Y (n x m)
*          int    n,m       I   size of matrix A,Y
*          double *X        O   X=A\Y or X=A'\Y (n x m)
* return : status (0:ok,0>:error)
* notes  : matirix stored by column-major order (fortran convention)
*          X can be same as Y
*-----------------------------------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info,*ipiv=imat(n,1);
    
    matcpy(B,A,n,n);
    matcpy(X,Y,n,m);
    dgetrf_(&n,&n,B,&n,ipiv,&info);
    if (!info) dgetrs_((char *)tr,&n,&m,B,&n,ipiv,X,&n,&info);
    free(ipiv); free(B); 
    return info;
}

#else /* without LAPACK/BLAS or MKL */

/* multiply matrix -----------------------------------------------------------*/
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C)
{
    double d;
    int i,j,x,f=tr[0]=='N'?(tr[1]=='N'?1:2):(tr[1]=='N'?3:4);
    
    for (i=0;i<n;i++) for (j=0;j<k;j++) {
        d=0.0;
        switch (f) {
            case 1: for (x=0;x<m;x++) d+=A[i+x*n]*B[x+j*m]; break;
            case 2: for (x=0;x<m;x++) d+=A[i+x*n]*B[j+x*k]; break;
            case 3: for (x=0;x<m;x++) d+=A[x+i*m]*B[x+j*m]; break;
            case 4: for (x=0;x<m;x++) d+=A[x+i*m]*B[j+x*k]; break;
        }
        if (beta==0.0) C[i+j*n]=alpha*d; else C[i+j*n]=alpha*d+beta*C[i+j*n];
    }
}
/* LU decomposition ----------------------------------------------------------*/
static int od_rtk_rtkcmn_ludcmp(double *A, int n, int *indx, double *d)
{
    double big,s,tmp,*vv=mat(n,1);
    int i,imax=0,j,k;
    
    *d=1.0;
    for (i=0;i<n;i++) {
        big=0.0; for (j=0;j<n;j++) if ((tmp=fabs(A[i+j*n]))>big) big=tmp;
        if (big>0.0) vv[i]=1.0/big; else {free(vv); return -1;}
    }
    for (j=0;j<n;j++) {
        for (i=0;i<j;i++) {
            s=A[i+j*n]; for (k=0;k<i;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
        }
        big=0.0;
        for (i=j;i<n;i++) {
            s=A[i+j*n]; for (k=0;k<j;k++) s-=A[i+k*n]*A[k+j*n]; A[i+j*n]=s;
            if ((tmp=vv[i]*fabs(s))>=big) {big=tmp; imax=i;}
        }
        if (j!=imax) {
            for (k=0;k<n;k++) {
                tmp=A[imax+k*n]; A[imax+k*n]=A[j+k*n]; A[j+k*n]=tmp;
            }
            *d=-(*d); vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (A[j+j*n]==0.0) {free(vv); return -1;}
        if (j!=n-1) {
            tmp=1.0/A[j+j*n]; for (i=j+1;i<n;i++) A[i+j*n]*=tmp;
        }
    }
    free(vv);
    return 0;
}
/* LU back-substitution ------------------------------------------------------*/
static void od_rtk_rtkcmn_lubksb(const double *A, int n, const int *indx, double *b)
{
    double s;
    int i,ii=-1,ip,j;
    
    for (i=0;i<n;i++) {
        ip=indx[i]; s=b[ip]; b[ip]=b[i];
        if (ii>=0) for (j=ii;j<i;j++) s-=A[i+j*n]*b[j]; else if (s) ii=i;
        b[i]=s;
    }
    for (i=n-1;i>=0;i--) {
        s=b[i]; for (j=i+1;j<n;j++) s-=A[i+j*n]*b[j]; b[i]=s/A[i+i*n];
    }
}
/* inverse of matrix ---------------------------------------------------------*/
extern int matinv(double *A, int n)
{
    double d,*B;
    int i,j,*indx;
    
    indx=imat(n,1); B=mat(n,n); matcpy(B,A,n,n);
    if (od_rtk_rtkcmn_ludcmp(B,n,indx,&d)) {free(indx); free(B); return -1;}
    for (j=0;j<n;j++) {
        for (i=0;i<n;i++) A[i+j*n]=0.0;
        A[j+j*n]=1.0;
        od_rtk_rtkcmn_lubksb(B,n,indx,A+j*n);
    }
    free(indx); free(B);
    return 0;
}
/* solve linear equation -----------------------------------------------------*/
extern int solve(const char *tr, const double *A, const double *Y, int n,
                 int m, double *X)
{
    double *B=mat(n,n);
    int info;
    
    matcpy(B,A,n,n);
    if (!(info=matinv(B,n))) matmul(tr[0]=='N'?"NN":"TN",n,m,n,1.0,B,Y,0.0,X);
    free(B);
    return info;
}
#endif
/* end of matrix routines ----------------------------------------------------*/

/* least square estimation -----------------------------------------------------
* least square estimation by solving normal equation (x=(A*A')^-1*A*y)
* args   : double *A        I   transpose of (weighted) design matrix (n x m)
*          double *y        I   (weighted) measurements (m x 1)
*          int    n,m       I   number of parameters and measurements (n<=m)
*          double *x        O   estmated parameters (n x 1)
*          double *Q        O   esimated parameters covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : for weighted least square, replace A and y by A*w and w*y (w=W^(1/2))
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q)
{
    double *Ay;
    int info;
    
    if (m<n) return -1;
    Ay=mat(n,1);
    matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
    matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
    if (!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}
/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
static int od_rtk_rtkcmn_filter_(const double *x, const double *P, const double *H,
                   const double *v, const double *R, int n, int m,
                   double *xp, double *Pp)
{
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    int info;
    
    matcpy(Q,R,m,m);
    matcpy(xp,x,n,1);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
    free(F); free(Q); free(K); free(I);
    return info;
}
extern int filter(double *x, double *P, const double *H, const double *v,
                  const double *R, int n, int m)
{
    double *x_,*xp_,*P_,*Pp_,*H_;
    int i,j,k,info,*ix;
    
    ix=imat(n,1); for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
    x_=mat(k,1); xp_=mat(k,1); P_=mat(k,k); Pp_=mat(k,k); H_=mat(k,m);
    for (i=0;i<k;i++) {
        x_[i]=x[ix[i]];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
    }
    info=od_rtk_rtkcmn_filter_(x_,P_,H_,v,R,k,m,xp_,Pp_);
    for (i=0;i<k;i++) {
        x[ix[i]]=xp_[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
    }
    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
    return info;
}
/* smoother --------------------------------------------------------------------
* combine forward and backward filters by fixed-interval smoother as follows:
*
*   xs=Qs*(Qf^-1*xf+Qb^-1*xb), Qs=(Qf^-1+Qb^-1)^-1)
*
* args   : double *xf       I   forward solutions (n x 1)
* args   : double *Qf       I   forward solutions covariance matrix (n x n)
*          double *xb       I   backward solutions (n x 1)
*          double *Qb       I   backward solutions covariance matrix (n x n)
*          int    n         I   number of solutions
*          double *xs       O   smoothed solutions (n x 1)
*          double *Qs       O   smoothed solutions covariance matrix (n x n)
* return : status (0:ok,0>:error)
* notes  : see reference [4] 5.2
*          matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern int smoother(const double *xf, const double *Qf, const double *xb,
                    const double *Qb, int n, double *xs, double *Qs)
{
    double *invQf=mat(n,n),*invQb=mat(n,n),*xx=mat(n,1);
    int i,info=-1;
    
    matcpy(invQf,Qf,n,n);
    matcpy(invQb,Qb,n,n);
    if (!matinv(invQf,n)&&!matinv(invQb,n)) {
        for (i=0;i<n*n;i++) Qs[i]=invQf[i]+invQb[i];
        if (!(info=matinv(Qs,n))) {
            matmul("NN",n,1,n,1.0,invQf,xf,0.0,xx);
            matmul("NN",n,1,n,1.0,invQb,xb,1.0,xx);
            matmul("NN",n,1,n,1.0,Qs,xx,0.0,xs);
        }
    }
    free(invQf); free(invQb); free(xx);
    return info;
}
/* print matrix ----------------------------------------------------------------
* print matrix to stdout
* args   : double *A        I   matrix A (n x m)
*          int    n,m       I   number of rows and columns of A
*          int    p,q       I   total columns, columns under decimal point
*         (FILE  *fp        I   output file pointer)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void matfprint(const double A[], int n, int m, int p, int q, FILE *fp)
{
    int i,j;
    
    for (i=0;i<n;i++) {
        for (j=0;j<m;j++) fprintf(fp," %*.*f",p,q,A[i+j*n]);
        fprintf(fp,"\n");
    }
}
extern void matprint(const double A[], int n, int m, int p, int q)
{
    matfprint(A,n,m,p,q,stdout);
}
/* string to number ------------------------------------------------------------
* convert substring in string to number
* args   : char   *s        I   string ("... nnn.nnn ...")
*          int    i,n       I   substring position and width
* return : converted number (0.0:error)
*-----------------------------------------------------------------------------*/
extern double str2num(const char *s, int i, int n)
{
    double value;
    char str[256],*p=str;
    
    if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<n) return 0.0;
    for (s+=i;*s&&--n>=0;s++) *p++=*s=='d'||*s=='D'?'E':*s;
    *p='\0';
    return sscanf(str,"%lf",&value)==1?value:0.0;
}
/* string to time --------------------------------------------------------------
* convert substring in string to gtime_t struct
* args   : char   *s        I   string ("... yyyy mm dd hh mm ss ...")
*          int    i,n       I   substring position and width
*          gtime_t *t       O   gtime_t struct
* return : status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int str2time(const char *s, int i, int n, gtime_t *t)
{
    double ep[6];
    char str[256],*p=str;
    
    if (i<0||(int)strlen(s)<i||(int)sizeof(str)-1<i) return -1;
    for (s+=i;*s&&--n>=0;) *p++=*s++;
    *p='\0';
    if (sscanf(str,"%lf %lf %lf %lf %lf %lf",ep,ep+1,ep+2,ep+3,ep+4,ep+5)<6)
        return -1;
    if (ep[0]<100.0) ep[0]+=ep[0]<80.0?2000.0:1900.0;
    *t=epoch2time(ep);
    return 0;
}
/* convert calendar day/time to time -------------------------------------------
* convert calendar day/time to gtime_t struct
* args   : double *ep       I   day/time {year,month,day,hour,min,sec}
* return : gtime_t struct
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern gtime_t epoch2time(const double *ep)
{
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    gtime_t time={0};
    int days,sec,year=(int)ep[0],mon=(int)ep[1],day=(int)ep[2];
    
    if (year<1970||2099<year||mon<1||12<mon) return time;
    
    /* leap year if year%4==0 in 1901-2099 */
    days=(year-1970)*365+(year-1969)/4+doy[mon-1]+day-2+(year%4==0&&mon>=3?1:0);
    sec=(int)floor(ep[5]);
    time.time=(time_t)days*86400+(int)ep[3]*3600+(int)ep[4]*60+sec;
    time.sec=ep[5]-sec;
    return time;
}
/* time to calendar day/time ---------------------------------------------------
* convert gtime_t struct to calendar day/time
* args   : gtime_t t        I   gtime_t struct
*          double *ep       O   day/time {year,month,day,hour,min,sec}
* return : none
* notes  : proper in 1970-2037 or 1970-2099 (64bit time_t)
*-----------------------------------------------------------------------------*/
extern void time2epoch(gtime_t t, double *ep)
{
    const int mday[]={ /* # of days in a month */
        31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
        31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
    };
    int days,sec,mon,day;
    
    /* leap year if year%4==0 in 1901-2099 */
    days=(int)(t.time/86400);
    sec=(int)(t.time-(time_t)days*86400);
    for (day=days%1461,mon=0;mon<48;mon++) {
        if (day>=mday[mon]) day-=mday[mon]; else break;
    }
    ep[0]=1970+days/1461*4+mon/12; ep[1]=mon%12+1; ep[2]=day+1;
    ep[3]=sec/3600; ep[4]=sec%3600/60; ep[5]=sec%60+t.sec;
}
/* gps time to time ------------------------------------------------------------
* convert week and tow in gps time to gtime_t struct
* args   : int    week      I   week number in gps time
*          double sec       I   time of week in gps time (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2time(int week, double sec)
{
    gtime_t t=epoch2time(od_rtk_rtkcmn_gpst0);
    
    if (sec<-1E9||1E9<sec) sec=0.0;
    t.time+=(time_t)86400*7*week+(int)sec;
    t.sec=sec-(int)sec;
    return t;
}
/* time to gps time ------------------------------------------------------------
* convert gtime_t struct to week and tow in gps time
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gps time (NULL: no output)
* return : time of week in gps time (s)
*-----------------------------------------------------------------------------*/
extern double time2gpst(gtime_t t, int *week)
{
    gtime_t t0=epoch2time(od_rtk_rtkcmn_gpst0);
    time_t sec=t.time-t0.time;
    int w=(int)(sec/(86400*7));
    
    if (week) *week=w;
    return (double)(sec-(double)w*86400*7)+t.sec;
}
/* galileo system time to time -------------------------------------------------
* convert week and tow in galileo system time (gst) to gtime_t struct
* args   : int    week      I   week number in gst
*          double sec       I   time of week in gst (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t gst2time(int week, double sec)
{
    gtime_t t=epoch2time(od_rtk_rtkcmn_gst0);
    
    if (sec<-1E9||1E9<sec) sec=0.0;
    t.time+=(time_t)86400*7*week+(int)sec;
    t.sec=sec-(int)sec;
    return t;
}
/* time to galileo system time -------------------------------------------------
* convert gtime_t struct to week and tow in galileo system time (gst)
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in gst (NULL: no output)
* return : time of week in gst (s)
*-----------------------------------------------------------------------------*/
extern double time2gst(gtime_t t, int *week)
{
    gtime_t t0=epoch2time(od_rtk_rtkcmn_gst0);
    time_t sec=t.time-t0.time;
    int w=(int)(sec/(86400*7));
    
    if (week) *week=w;
    return (double)(sec-(double)w*86400*7)+t.sec;
}
/* beidou time (bdt) to time ---------------------------------------------------
* convert week and tow in beidou time (bdt) to gtime_t struct
* args   : int    week      I   week number in bdt
*          double sec       I   time of week in bdt (s)
* return : gtime_t struct
*-----------------------------------------------------------------------------*/
extern gtime_t bdt2time(int week, double sec)
{
    gtime_t t=epoch2time(od_rtk_rtkcmn_bdt0);
    
    if (sec<-1E9||1E9<sec) sec=0.0;
    t.time+=(time_t)86400*7*week+(int)sec;
    t.sec=sec-(int)sec;
    return t;
}
/* time to beidouo time (bdt) --------------------------------------------------
* convert gtime_t struct to week and tow in beidou time (bdt)
* args   : gtime_t t        I   gtime_t struct
*          int    *week     IO  week number in bdt (NULL: no output)
* return : time of week in bdt (s)
*-----------------------------------------------------------------------------*/
extern double time2bdt(gtime_t t, int *week)
{
    gtime_t t0=epoch2time(od_rtk_rtkcmn_bdt0);
    time_t sec=t.time-t0.time;
    int w=(int)(sec/(86400*7));
    
    if (week) *week=w;
    return (double)(sec-(double)w*86400*7)+t.sec;
}
/* add time --------------------------------------------------------------------
* add time to gtime_t struct
* args   : gtime_t t        I   gtime_t struct
*          double sec       I   time to add (s)
* return : gtime_t struct (t+sec)
*-----------------------------------------------------------------------------*/
extern gtime_t timeadd(gtime_t t, double sec)
{
    double tt;
    
    t.sec+=sec; tt=floor(t.sec); t.time+=(int)tt; t.sec-=tt;
    return t;
}
/* time difference -------------------------------------------------------------
* difference between gtime_t structs
* args   : gtime_t t1,t2    I   gtime_t structs
* return : time difference (t1-t2) (s)
*-----------------------------------------------------------------------------*/
extern double timediff(gtime_t t1, gtime_t t2)
{
    return difftime(t1.time,t2.time)+t1.sec-t2.sec;
}
/* get current time in utc -----------------------------------------------------
* get current time in utc
* args   : none
* return : current time in utc
*-----------------------------------------------------------------------------*/
static double od_rtk_rtkcmn_timeoffset_=0.0;        /* time offset (s) */

extern gtime_t timeget(void)
{
    gtime_t time;
    double ep[6]={0};
#ifdef WIN32
    SYSTEMTIME ts;
    
    GetSystemTime(&ts); /* utc */
    ep[0]=ts.wYear; ep[1]=ts.wMonth;  ep[2]=ts.wDay;
    ep[3]=ts.wHour; ep[4]=ts.wMinute; ep[5]=ts.wSecond+ts.wMilliseconds*1E-3;
#else
    struct timeval tv;
    struct tm *tt;
    
    if (!gettimeofday(&tv,NULL)&&(tt=gmtime(&tv.tv_sec))) {
        ep[0]=tt->tm_year+1900; ep[1]=tt->tm_mon+1; ep[2]=tt->tm_mday;
        ep[3]=tt->tm_hour; ep[4]=tt->tm_min; ep[5]=tt->tm_sec+tv.tv_usec*1E-6;
    }
#endif
    time=epoch2time(ep);
    
#ifdef CPUTIME_IN_GPST /* cputime operated in gpst */
    time=gpst2utc(time);
#endif
    return timeadd(time,od_rtk_rtkcmn_timeoffset_);
}
/* set current time in utc -----------------------------------------------------
* set current time in utc
* args   : gtime_t          I   current time in utc
* return : none
* notes  : just set time offset between cpu time and current time
*          the time offset is reflected to only timeget()
*          not reentrant
*-----------------------------------------------------------------------------*/
extern void timeset(gtime_t t)
{
    od_rtk_rtkcmn_timeoffset_+=timediff(t,timeget());
}
/* reset current time ----------------------------------------------------------
* reset current time
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void timereset(void)
{
    od_rtk_rtkcmn_timeoffset_=0.0;
}
/* read leap seconds table by text -------------------------------------------*/
static int od_rtk_rtkcmn_read_leaps_text(FILE *fp)
{
    char buff[256],*p;
    int i,n=0,ep[6],ls;
    
    rewind(fp);
    
    while (fgets(buff,sizeof(buff),fp)&&n<MAXLEAPS) {
        if ((p=strchr(buff,'#'))) *p='\0';
        if (sscanf(buff,"%d %d %d %d %d %d %d",ep,ep+1,ep+2,ep+3,ep+4,ep+5,
                   &ls)<7) continue;
        for (i=0;i<6;i++) od_rtk_rtkcmn_leaps[n][i]=ep[i];
        od_rtk_rtkcmn_leaps[n++][6]=ls;
    }
    return n;
}
/* read leap seconds table by usno -------------------------------------------*/
static int od_rtk_rtkcmn_read_leaps_usno(FILE *fp)
{
    static const char *months[]={
        "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"
    };
    double jd,tai_utc;
    char buff[256],month[32],ls[MAXLEAPS][7]={{0}};
    int i,j,y,m,d,n=0;
    
    rewind(fp);
    
    while (fgets(buff,sizeof(buff),fp)&&n<MAXLEAPS) {
        if (sscanf(buff,"%d %s %d =JD %lf TAI-UTC= %lf",&y,month,&d,&jd,
                   &tai_utc)<5) continue;
        if (y<1980) continue;
        for (m=1;m<=12;m++) if (!strcmp(months[m-1],month)) break;
        if (m>=13) continue;
        ls[n][0]=y;
        ls[n][1]=m;
        ls[n][2]=d;
        ls[n++][6]=(char)(19.0-tai_utc);
    }
    for (i=0;i<n;i++) for (j=0;j<7;j++) {
        od_rtk_rtkcmn_leaps[i][j]=ls[n-i-1][j];
    }
    return n;
}
/* read leap seconds table -----------------------------------------------------
* read leap seconds table
* args   : char    *file    I   leap seconds table file
* return : status (1:ok,0:error)
* notes  : The leap second table should be as follows or leapsec.dat provided
*          by USNO.
*          (1) The records in the table file cosist of the following fields:
*              year month day hour min sec UTC-GPST(s)
*          (2) The date and time indicate the start UTC time for the UTC-GPST
*          (3) The date and time should be descending order.
*-----------------------------------------------------------------------------*/
extern int read_leaps(const char *file)
{
    FILE *fp;
    int i,n;
    
    if (!(fp=fopen(file,"r"))) return 0;
    
    /* read leap seconds table by text or usno */
    if (!(n=od_rtk_rtkcmn_read_leaps_text(fp))&&!(n=od_rtk_rtkcmn_read_leaps_usno(fp))) {
        fclose(fp);
        return 0;
    }
    for (i=0;i<7;i++) od_rtk_rtkcmn_leaps[n][i]=0.0;
    fclose(fp);
    return 1;
}
/* gpstime to utc --------------------------------------------------------------
* convert gpstime to utc considering leap seconds
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in utc
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2utc(gtime_t t)
{
    gtime_t tu;
    int i;
    
    for (i=0;od_rtk_rtkcmn_leaps[i][0]>0;i++) {
        tu=timeadd(t,od_rtk_rtkcmn_leaps[i][6]);
        if (timediff(tu,epoch2time(od_rtk_rtkcmn_leaps[i]))>=0.0) return tu;
    }
    return t;
}
/* utc to gpstime --------------------------------------------------------------
* convert utc to gpstime considering leap seconds
* args   : gtime_t t        I   time expressed in utc
* return : time expressed in gpstime
* notes  : ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
extern gtime_t utc2gpst(gtime_t t)
{
    int i;
    
    for (i=0;od_rtk_rtkcmn_leaps[i][0]>0;i++) {
        if (timediff(t,epoch2time(od_rtk_rtkcmn_leaps[i]))>=0.0) return timeadd(t,-od_rtk_rtkcmn_leaps[i][6]);
    }
    return t;
}
/* gpstime to bdt --------------------------------------------------------------
* convert gpstime to bdt (beidou navigation satellite system time)
* args   : gtime_t t        I   time expressed in gpstime
* return : time expressed in bdt
* notes  : ref [8] 3.3, 2006/1/1 00:00 BDT = 2006/1/1 00:00 UTC
*          no leap seconds in BDT
*          ignore slight time offset under 100 ns
*-----------------------------------------------------------------------------*/
extern gtime_t gpst2bdt(gtime_t t)
{
    return timeadd(t,-14.0);
}
/* bdt to gpstime --------------------------------------------------------------
* convert bdt (beidou navigation satellite system time) to gpstime
* args   : gtime_t t        I   time expressed in bdt
* return : time expressed in gpstime
* notes  : see gpst2bdt()
*-----------------------------------------------------------------------------*/
extern gtime_t bdt2gpst(gtime_t t)
{
    return timeadd(t,14.0);
}
/* time to day and sec -------------------------------------------------------*/
static double od_rtk_rtkcmn_time2sec(gtime_t time, gtime_t *day)
{
    double ep[6],sec;
    time2epoch(time,ep);
    sec=ep[3]*3600.0+ep[4]*60.0+ep[5];
    ep[3]=ep[4]=ep[5]=0.0;
    *day=epoch2time(ep);
    return sec;
}
/* utc to gmst -----------------------------------------------------------------
* convert utc to gmst (Greenwich mean sidereal time)
* args   : gtime_t t        I   time expressed in utc
*          double ut1_utc   I   UT1-UTC (s)
* return : gmst (rad)
*-----------------------------------------------------------------------------*/
extern double utc2gmst(gtime_t t, double ut1_utc)
{
    const double ep2000[]={2000,1,1,12,0,0};
    gtime_t tut,tut0;
    double ut,t1,t2,t3,gmst0,gmst;
    
    tut=timeadd(t,ut1_utc);
    ut=od_rtk_rtkcmn_time2sec(tut,&tut0);
    t1=timediff(tut0,epoch2time(ep2000))/86400.0/36525.0;
    t2=t1*t1; t3=t2*t1;
    gmst0=24110.54841+8640184.812866*t1+0.093104*t2-6.2E-6*t3;
    gmst=gmst0+1.002737909350795*ut;
    
    return fmod(gmst,86400.0)*PI/43200.0; /* 0 <= gmst <= 2*PI */
}
/* time to string --------------------------------------------------------------
* convert gtime_t struct to string
* args   : gtime_t t        I   gtime_t struct
*          char   *s        O   string ("yyyy/mm/dd hh:mm:ss.ssss")
*          int    n         I   number of decimals
* return : none
*-----------------------------------------------------------------------------*/
extern void time2str(gtime_t t, char *s, int n)
{
    double ep[6];
    
    if (n<0) n=0; else if (n>12) n=12;
    if (1.0-t.sec<0.5/pow(10.0,n)) {t.time++; t.sec=0.0;};
    time2epoch(t,ep);
    sprintf(s,"%04.0f/%02.0f/%02.0f %02.0f:%02.0f:%0*.*f",ep[0],ep[1],ep[2],
            ep[3],ep[4],n<=0?2:n+3,n<=0?0:n,ep[5]);
}
/* get time string -------------------------------------------------------------
* get time string
* args   : gtime_t t        I   gtime_t struct
*          int    n         I   number of decimals
* return : time string
* notes  : not reentrant, do not use multiple in a function
*-----------------------------------------------------------------------------*/
extern char *time_str(gtime_t t, int n)
{
    static char buff[64];
    time2str(t,buff,n);
    return buff;
}
/* time to day of year ---------------------------------------------------------
* convert time to day of year
* args   : gtime_t t        I   gtime_t struct
* return : day of year (days)
*-----------------------------------------------------------------------------*/
extern double time2doy(gtime_t t)
{
    double ep[6];
    
    time2epoch(t,ep);
    ep[1]=ep[2]=1.0; ep[3]=ep[4]=ep[5]=0.0;
    return timediff(t,epoch2time(ep))/86400.0+1.0;
}
/* adjust gps week number ------------------------------------------------------
* adjust gps week number using cpu time
* args   : int   week       I   not-adjusted gps week number (0-1023)
* return : adjusted gps week number
*-----------------------------------------------------------------------------*/
extern int adjgpsweek(int week)
{
    int w;
    (void)time2gpst(utc2gpst(timeget()),&w);
    if (w<1560) w=1560; /* use 2009/12/1 if time is earlier than 2009/12/1 */
    return week+(w-week+1)/1024*1024;
}
/* get tick time ---------------------------------------------------------------
* get current tick in ms
* args   : none
* return : current tick in ms
*-----------------------------------------------------------------------------*/
extern uint32_t tickget(void)
{
#ifdef WIN32
    return (uint32_t)timeGetTime();
#else
    struct timespec tp={0};
    struct timeval  tv={0};
    
#ifdef CLOCK_MONOTONIC_RAW
    /* linux kernel > 2.6.28 */
    if (!clock_gettime(CLOCK_MONOTONIC_RAW,&tp)) {
        return tp.tv_sec*1000u+tp.tv_nsec/1000000u;
    }
    else {
        gettimeofday(&tv,NULL);
        return tv.tv_sec*1000u+tv.tv_usec/1000u;
    }
#else
    gettimeofday(&tv,NULL);
    return tv.tv_sec*1000u+tv.tv_usec/1000u;
#endif
#endif /* WIN32 */
}
/* sleep ms --------------------------------------------------------------------
* sleep ms
* args   : int   ms         I   miliseconds to sleep (<0:no sleep)
* return : none
*-----------------------------------------------------------------------------*/
extern void sleepms(int ms)
{
#ifdef WIN32
    if (ms<5) Sleep(1); else Sleep(ms);
#else
    struct timespec ts;
    if (ms<=0) return;
    ts.tv_sec=(time_t)(ms/1000);
    ts.tv_nsec=(long)(ms%1000*1000000);
    nanosleep(&ts,NULL);
#endif
}
/* convert degree to deg-min-sec -----------------------------------------------
* convert degree to degree-minute-second
* args   : double deg       I   degree
*          double *dms      O   degree-minute-second {deg,min,sec}
*          int    ndec      I   number of decimals of second
* return : none
*-----------------------------------------------------------------------------*/
extern void deg2dms(double deg, double *dms, int ndec)
{
    double sign=deg<0.0?-1.0:1.0,a=fabs(deg);
    double unit=pow(0.1,ndec);
    dms[0]=floor(a); a=(a-dms[0])*60.0;
    dms[1]=floor(a); a=(a-dms[1])*60.0;
    dms[2]=floor(a/unit+0.5)*unit;
    if (dms[2]>=60.0) {
        dms[2]=0.0;
        dms[1]+=1.0;
        if (dms[1]>=60.0) {
            dms[1]=0.0;
            dms[0]+=1.0;
        }
    }
    dms[0]*=sign;
}
/* convert deg-min-sec to degree -----------------------------------------------
* convert degree-minute-second to degree
* args   : double *dms      I   degree-minute-second {deg,min,sec}
* return : degree
*-----------------------------------------------------------------------------*/
extern double dms2deg(const double *dms)
{
    double sign=dms[0]<0.0?-1.0:1.0;
    return sign*(fabs(dms[0])+dms[1]/60.0+dms[2]/3600.0);
}
/* transform ecef to geodetic postion ------------------------------------------
* transform ecef position to geodetic position
* args   : double *r        I   ecef position {x,y,z} (m)
*          double *pos      O   geodetic position {lat,lon,h} (rad,m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos)
{
    double e2=FE_WGS84*(2.0-FE_WGS84),r2=dot(r,r,2),z,zk,v=RE_WGS84,sinp;
    
    for (z=r[2],zk=0.0;fabs(z-zk)>=1E-4;) {
        zk=z;
        sinp=z/sqrt(r2+z*z);
        v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
        z=r[2]+v*e2*sinp;
    }
    pos[0]=r2>1E-12?atan(z/sqrt(r2)):(r[2]>0.0?PI/2.0:-PI/2.0);
    pos[1]=r2>1E-12?atan2(r[1],r[0]):0.0;
    pos[2]=sqrt(r2+z*z)-v;
}
/* transform geodetic to ecef position -----------------------------------------
* transform geodetic position to ecef position
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *r        O   ecef position {x,y,z} (m)
* return : none
* notes  : WGS84, ellipsoidal height
*-----------------------------------------------------------------------------*/
extern void pos2ecef(const double *pos, double *r)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    double e2=FE_WGS84*(2.0-FE_WGS84),v=RE_WGS84/sqrt(1.0-e2*sinp*sinp);
    
    r[0]=(v+pos[2])*cosp*cosl;
    r[1]=(v+pos[2])*cosp*sinl;
    r[2]=(v*(1.0-e2)+pos[2])*sinp;
}
/* ecef to local coordinate transfromation matrix ------------------------------
* compute ecef to local coordinate transfromation matrix
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *E        O   ecef to local coord transformation matrix (3x3)
* return : none
* notes  : matirix stored by column-major order (fortran convention)
*-----------------------------------------------------------------------------*/
extern void xyz2enu(const double *pos, double *E)
{
    double sinp=sin(pos[0]),cosp=cos(pos[0]),sinl=sin(pos[1]),cosl=cos(pos[1]);
    
    E[0]=-sinl;      E[3]=cosl;       E[6]=0.0;
    E[1]=-sinp*cosl; E[4]=-sinp*sinl; E[7]=cosp;
    E[2]=cosp*cosl;  E[5]=cosp*sinl;  E[8]=sinp;
}
/* transform ecef vector to local tangental coordinate -------------------------
* transform ecef vector to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *r        I   vector in ecef coordinate {x,y,z}
*          double *e        O   vector in local tangental coordinate {e,n,u}
* return : none
*-----------------------------------------------------------------------------*/
extern void ecef2enu(const double *pos, const double *r, double *e)
{
    double E[9];
    
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,r,0.0,e);
}
/* transform local vector to ecef coordinate -----------------------------------
* transform local tangental coordinate vector to ecef
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *e        I   vector in local tangental coordinate {e,n,u}
*          double *r        O   vector in ecef coordinate {x,y,z}
* return : none
*-----------------------------------------------------------------------------*/
extern void enu2ecef(const double *pos, const double *e, double *r)
{
    double E[9];
    
    xyz2enu(pos,E);
    matmul("TN",3,1,3,1.0,E,e,0.0,r);
}
/* transform covariance to local tangental coordinate --------------------------
* transform ecef covariance to local tangental coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *P        I   covariance in ecef coordinate
*          double *Q        O   covariance in local tangental coordinate
* return : none
*-----------------------------------------------------------------------------*/
extern void covenu(const double *pos, const double *P, double *Q)
{
    double E[9],EP[9];
    
    xyz2enu(pos,E);
    matmul("NN",3,3,3,1.0,E,P,0.0,EP);
    matmul("NT",3,3,3,1.0,EP,E,0.0,Q);
}
/* transform local enu coordinate covariance to xyz-ecef -----------------------
* transform local enu covariance to xyz-ecef coordinate
* args   : double *pos      I   geodetic position {lat,lon} (rad)
*          double *Q        I   covariance in local enu coordinate
*          double *P        O   covariance in xyz-ecef coordinate
* return : none
*-----------------------------------------------------------------------------*/
extern void covecef(const double *pos, const double *Q, double *P)
{
    double E[9],EQ[9];
    
    xyz2enu(pos,E);
    matmul("TN",3,3,3,1.0,E,Q,0.0,EQ);
    matmul("NN",3,3,3,1.0,EQ,E,0.0,P);
}
/* coordinate rotation matrix ------------------------------------------------*/
#define Rx(t,X) do { \
    (X)[0]=1.0; (X)[1]=(X)[2]=(X)[3]=(X)[6]=0.0; \
    (X)[4]=(X)[8]=cos(t); (X)[7]=sin(t); (X)[5]=-(X)[7]; \
} while (0)

#define Ry(t,X) do { \
    (X)[4]=1.0; (X)[1]=(X)[3]=(X)[5]=(X)[7]=0.0; \
    (X)[0]=(X)[8]=cos(t); (X)[2]=sin(t); (X)[6]=-(X)[2]; \
} while (0)

#define Rz(t,X) do { \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

/* astronomical arguments: f={l,l',F,D,OMG} (rad) ----------------------------*/
static void od_rtk_rtkcmn_ast_args(double t, double *f)
{
    static const double fc[][5]={ /* coefficients for iau 1980 nutation */
        { 134.96340251, 1717915923.2178,  31.8792,  0.051635, -0.00024470},
        { 357.52910918,  129596581.0481,  -0.5532,  0.000136, -0.00001149},
        {  93.27209062, 1739527262.8478, -12.7512, -0.001037,  0.00000417},
        { 297.85019547, 1602961601.2090,  -6.3706,  0.006593, -0.00003169},
        { 125.04455501,   -6962890.2665,   7.4722,  0.007702, -0.00005939}
    };
    double tt[4];
    int i,j;
    
    for (tt[0]=t,i=1;i<4;i++) tt[i]=tt[i-1]*t;
    for (i=0;i<5;i++) {
        f[i]=fc[i][0]*3600.0;
        for (j=0;j<4;j++) f[i]+=fc[i][j+1]*tt[j];
        f[i]=fmod(f[i]*AS2R,2.0*PI);
    }
}
/* iau 1980 nutation ---------------------------------------------------------*/
static void od_rtk_rtkcmn_nut_iau1980(double t, const double *f, double *dpsi, double *deps)
{
    static const double nut[106][10]={
        {   0,   0,   0,   0,   1, -6798.4, -171996, -174.2, 92025,   8.9},
        {   0,   0,   2,  -2,   2,   182.6,  -13187,   -1.6,  5736,  -3.1},
        {   0,   0,   2,   0,   2,    13.7,   -2274,   -0.2,   977,  -0.5},
        {   0,   0,   0,   0,   2, -3399.2,    2062,    0.2,  -895,   0.5},
        {   0,  -1,   0,   0,   0,  -365.3,   -1426,    3.4,    54,  -0.1},
        {   1,   0,   0,   0,   0,    27.6,     712,    0.1,    -7,   0.0},
        {   0,   1,   2,  -2,   2,   121.7,    -517,    1.2,   224,  -0.6},
        {   0,   0,   2,   0,   1,    13.6,    -386,   -0.4,   200,   0.0},
        {   1,   0,   2,   0,   2,     9.1,    -301,    0.0,   129,  -0.1},
        {   0,  -1,   2,  -2,   2,   365.2,     217,   -0.5,   -95,   0.3},
        {  -1,   0,   0,   2,   0,    31.8,     158,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   1,   177.8,     129,    0.1,   -70,   0.0},
        {  -1,   0,   2,   0,   2,    27.1,     123,    0.0,   -53,   0.0},
        {   1,   0,   0,   0,   1,    27.7,      63,    0.1,   -33,   0.0},
        {   0,   0,   0,   2,   0,    14.8,      63,    0.0,    -2,   0.0},
        {  -1,   0,   2,   2,   2,     9.6,     -59,    0.0,    26,   0.0},
        {  -1,   0,   0,   0,   1,   -27.4,     -58,   -0.1,    32,   0.0},
        {   1,   0,   2,   0,   1,     9.1,     -51,    0.0,    27,   0.0},
        {  -2,   0,   0,   2,   0,  -205.9,     -48,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   1,  1305.5,      46,    0.0,   -24,   0.0},
        {   0,   0,   2,   2,   2,     7.1,     -38,    0.0,    16,   0.0},
        {   2,   0,   2,   0,   2,     6.9,     -31,    0.0,    13,   0.0},
        {   2,   0,   0,   0,   0,    13.8,      29,    0.0,    -1,   0.0},
        {   1,   0,   2,  -2,   2,    23.9,      29,    0.0,   -12,   0.0},
        {   0,   0,   2,   0,   0,    13.6,      26,    0.0,    -1,   0.0},
        {   0,   0,   2,  -2,   0,   173.3,     -22,    0.0,     0,   0.0},
        {  -1,   0,   2,   0,   1,    27.0,      21,    0.0,   -10,   0.0},
        {   0,   2,   0,   0,   0,   182.6,      17,   -0.1,     0,   0.0},
        {   0,   2,   2,  -2,   2,    91.3,     -16,    0.1,     7,   0.0},
        {  -1,   0,   0,   2,   1,    32.0,      16,    0.0,    -8,   0.0},
        {   0,   1,   0,   0,   1,   386.0,     -15,    0.0,     9,   0.0},
        {   1,   0,   0,  -2,   1,   -31.7,     -13,    0.0,     7,   0.0},
        {   0,  -1,   0,   0,   1,  -346.6,     -12,    0.0,     6,   0.0},
        {   2,   0,  -2,   0,   0, -1095.2,      11,    0.0,     0,   0.0},
        {  -1,   0,   2,   2,   1,     9.5,     -10,    0.0,     5,   0.0},
        {   1,   0,   2,   2,   2,     5.6,      -8,    0.0,     3,   0.0},
        {   0,  -1,   2,   0,   2,    14.2,      -7,    0.0,     3,   0.0},
        {   0,   0,   2,   2,   1,     7.1,      -7,    0.0,     3,   0.0},
        {   1,   1,   0,  -2,   0,   -34.8,      -7,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   2,    13.2,       7,    0.0,    -3,   0.0},
        {  -2,   0,   0,   2,   1,  -199.8,      -6,    0.0,     3,   0.0},
        {   0,   0,   0,   2,   1,    14.8,      -6,    0.0,     3,   0.0},
        {   2,   0,   2,  -2,   2,    12.8,       6,    0.0,    -3,   0.0},
        {   1,   0,   0,   2,   0,     9.6,       6,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   1,    23.9,       6,    0.0,    -3,   0.0},
        {   0,   0,   0,  -2,   1,   -14.7,      -5,    0.0,     3,   0.0},
        {   0,  -1,   2,  -2,   1,   346.6,      -5,    0.0,     3,   0.0},
        {   2,   0,   2,   0,   1,     6.9,      -5,    0.0,     3,   0.0},
        {   1,  -1,   0,   0,   0,    29.8,       5,    0.0,     0,   0.0},
        {   1,   0,   0,  -1,   0,   411.8,      -4,    0.0,     0,   0.0},
        {   0,   0,   0,   1,   0,    29.5,      -4,    0.0,     0,   0.0},
        {   0,   1,   0,  -2,   0,   -15.4,      -4,    0.0,     0,   0.0},
        {   1,   0,  -2,   0,   0,   -26.9,       4,    0.0,     0,   0.0},
        {   2,   0,   0,  -2,   1,   212.3,       4,    0.0,    -2,   0.0},
        {   0,   1,   2,  -2,   1,   119.6,       4,    0.0,    -2,   0.0},
        {   1,   1,   0,   0,   0,    25.6,      -3,    0.0,     0,   0.0},
        {   1,  -1,   0,  -1,   0, -3232.9,      -3,    0.0,     0,   0.0},
        {  -1,  -1,   2,   2,   2,     9.8,      -3,    0.0,     1,   0.0},
        {   0,  -1,   2,   2,   2,     7.2,      -3,    0.0,     1,   0.0},
        {   1,  -1,   2,   0,   2,     9.4,      -3,    0.0,     1,   0.0},
        {   3,   0,   2,   0,   2,     5.5,      -3,    0.0,     1,   0.0},
        {  -2,   0,   2,   0,   2,  1615.7,      -3,    0.0,     1,   0.0},
        {   1,   0,   2,   0,   0,     9.1,       3,    0.0,     0,   0.0},
        {  -1,   0,   2,   4,   2,     5.8,      -2,    0.0,     1,   0.0},
        {   1,   0,   0,   0,   2,    27.8,      -2,    0.0,     1,   0.0},
        {  -1,   0,   2,  -2,   1,   -32.6,      -2,    0.0,     1,   0.0},
        {   0,  -2,   2,  -2,   1,  6786.3,      -2,    0.0,     1,   0.0},
        {  -2,   0,   0,   0,   1,   -13.7,      -2,    0.0,     1,   0.0},
        {   2,   0,   0,   0,   1,    13.8,       2,    0.0,    -1,   0.0},
        {   3,   0,   0,   0,   0,     9.2,       2,    0.0,     0,   0.0},
        {   1,   1,   2,   0,   2,     8.9,       2,    0.0,    -1,   0.0},
        {   0,   0,   2,   1,   2,     9.3,       2,    0.0,    -1,   0.0},
        {   1,   0,   0,   2,   1,     9.6,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,   2,   1,     5.6,      -1,    0.0,     1,   0.0},
        {   1,   1,   0,  -2,   1,   -34.7,      -1,    0.0,     0,   0.0},
        {   0,   1,   0,   2,   0,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   1,   2,  -2,   0,   117.5,      -1,    0.0,     0,   0.0},
        {   0,   1,  -2,   2,   0,  -329.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,   2,   0,    23.8,      -1,    0.0,     0,   0.0},
        {   1,   0,  -2,  -2,   0,    -9.5,      -1,    0.0,     0,   0.0},
        {   1,   0,   2,  -2,   0,    32.8,      -1,    0.0,     0,   0.0},
        {   1,   0,   0,  -4,   0,   -10.1,      -1,    0.0,     0,   0.0},
        {   2,   0,   0,  -4,   0,   -15.9,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,   4,   2,     4.8,      -1,    0.0,     0,   0.0},
        {   0,   0,   2,  -1,   2,    25.4,      -1,    0.0,     0,   0.0},
        {  -2,   0,   2,   4,   2,     7.3,      -1,    0.0,     1,   0.0},
        {   2,   0,   2,   2,   2,     4.7,      -1,    0.0,     0,   0.0},
        {   0,  -1,   2,   0,   1,    14.2,      -1,    0.0,     0,   0.0},
        {   0,   0,  -2,   0,   1,   -13.6,      -1,    0.0,     0,   0.0},
        {   0,   0,   4,  -2,   2,    12.7,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   0,   2,   409.2,       1,    0.0,     0,   0.0},
        {   1,   1,   2,  -2,   2,    22.5,       1,    0.0,    -1,   0.0},
        {   3,   0,   2,  -2,   2,     8.7,       1,    0.0,     0,   0.0},
        {  -2,   0,   2,   2,   2,    14.6,       1,    0.0,    -1,   0.0},
        {  -1,   0,   0,   0,   2,   -27.3,       1,    0.0,    -1,   0.0},
        {   0,   0,  -2,   2,   1,  -169.0,       1,    0.0,     0,   0.0},
        {   0,   1,   2,   0,   1,    13.1,       1,    0.0,     0,   0.0},
        {  -1,   0,   4,   0,   2,     9.1,       1,    0.0,     0,   0.0},
        {   2,   1,   0,  -2,   0,   131.7,       1,    0.0,     0,   0.0},
        {   2,   0,   0,   2,   0,     7.1,       1,    0.0,     0,   0.0},
        {   2,   0,   2,  -2,   1,    12.8,       1,    0.0,    -1,   0.0},
        {   2,   0,  -2,   0,   1,  -943.2,       1,    0.0,     0,   0.0},
        {   1,  -1,   0,  -2,   0,   -29.3,       1,    0.0,     0,   0.0},
        {  -1,   0,   0,   1,   1,  -388.3,       1,    0.0,     0,   0.0},
        {  -1,  -1,   0,   2,   1,    35.0,       1,    0.0,     0,   0.0},
        {   0,   1,   0,   1,   0,    27.3,       1,    0.0,     0,   0.0}
    };
    double ang;
    int i,j;
    
    *dpsi=*deps=0.0;
    
    for (i=0;i<106;i++) {
        ang=0.0;
        for (j=0;j<5;j++) ang+=nut[i][j]*f[j];
        *dpsi+=(nut[i][6]+nut[i][7]*t)*sin(ang);
        *deps+=(nut[i][8]+nut[i][9]*t)*cos(ang);
    }
    *dpsi*=1E-4*AS2R; /* 0.1 mas -> rad */
    *deps*=1E-4*AS2R;
}
/* eci to ecef transformation matrix -------------------------------------------
* compute eci to ecef transformation matrix
* args   : gtime_t tutc     I   time in utc
*          double *erpv     I   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*          double *U        O   eci to ecef transformation matrix (3 x 3)
*          double *gmst     IO  greenwich mean sidereal time (rad)
*                               (NULL: no output)
* return : none
* note   : see ref [3] chap 5
*          not thread-safe
*-----------------------------------------------------------------------------*/
extern void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst)
{
    const double ep2000[]={2000,1,1,12,0,0};
    static gtime_t tutc_;
    static double U_[9],gmst_;
    gtime_t tgps;
    double eps,ze,th,z,t,t2,t3,dpsi,deps,gast,f[5];
    double R1[9],R2[9],R3[9],R[9],W[9],N[9],P[9],NP[9];
    int i;
    
    trace(4,"eci2ecef: tutc=%s\n",time_str(tutc,3));
    
    if (fabs(timediff(tutc,tutc_))<0.01) { /* read cache */
        for (i=0;i<9;i++) U[i]=U_[i];
        if (gmst) *gmst=gmst_; 
        return;
    }
    tutc_=tutc;
    
    /* terrestrial time */
    tgps=utc2gpst(tutc_);
    t=(timediff(tgps,epoch2time(ep2000))+19.0+32.184)/86400.0/36525.0;
    t2=t*t; t3=t2*t;
    
    /* astronomical arguments */
    od_rtk_rtkcmn_ast_args(t,f);
    
    /* iau 1976 precession */
    ze=(2306.2181*t+0.30188*t2+0.017998*t3)*AS2R;
    th=(2004.3109*t-0.42665*t2-0.041833*t3)*AS2R;
    z =(2306.2181*t+1.09468*t2+0.018203*t3)*AS2R;
    eps=(84381.448-46.8150*t-0.00059*t2+0.001813*t3)*AS2R;
    Rz(-z,R1); Ry(th,R2); Rz(-ze,R3);
    matmul("NN",3,3,3,1.0,R1,R2,0.0,R);
    matmul("NN",3,3,3,1.0,R, R3,0.0,P); /* P=Rz(-z)*Ry(th)*Rz(-ze) */
    
    /* iau 1980 nutation */
    od_rtk_rtkcmn_nut_iau1980(t,f,&dpsi,&deps);
    Rx(-eps-deps,R1); Rz(-dpsi,R2); Rx(eps,R3);
    matmul("NN",3,3,3,1.0,R1,R2,0.0,R);
    matmul("NN",3,3,3,1.0,R ,R3,0.0,N); /* N=Rx(-eps)*Rz(-dspi)*Rx(eps) */
    
    /* greenwich aparent sidereal time (rad) */
    gmst_=utc2gmst(tutc_,erpv[2]);
    gast=gmst_+dpsi*cos(eps);
    gast+=(0.00264*sin(f[4])+0.000063*sin(2.0*f[4]))*AS2R;
    
    /* eci to ecef transformation matrix */
    Ry(-erpv[0],R1); Rx(-erpv[1],R2); Rz(gast,R3);
    matmul("NN",3,3,3,1.0,R1,R2,0.0,W );
    matmul("NN",3,3,3,1.0,W ,R3,0.0,R ); /* W=Ry(-xp)*Rx(-yp) */
    matmul("NN",3,3,3,1.0,N ,P ,0.0,NP);
    matmul("NN",3,3,3,1.0,R ,NP,0.0,U_); /* U=W*Rz(gast)*N*P */
    
    for (i=0;i<9;i++) U[i]=U_[i];
    if (gmst) *gmst=gmst_; 
    
    trace(5,"gmst=%.12f gast=%.12f\n",gmst_,gast);
    trace(5,"P=\n"); tracemat(5,P,3,3,15,12);
    trace(5,"N=\n"); tracemat(5,N,3,3,15,12);
    trace(5,"W=\n"); tracemat(5,W,3,3,15,12);
    trace(5,"U=\n"); tracemat(5,U,3,3,15,12);
}
/* decode antenna parameter field --------------------------------------------*/
static int od_rtk_rtkcmn_decodef(char *p, int n, double *v)
{
    int i;
    
    for (i=0;i<n;i++) v[i]=0.0;
    for (i=0,p=strtok(p," ");p&&i<n;p=strtok(NULL," ")) {
        v[i++]=atof(p)*1E-3;
    }
    return i;
}
/* add antenna parameter -----------------------------------------------------*/
static void od_rtk_rtkcmn_addpcv(const pcv_t *pcv, pcvs_t *pcvs)
{
    pcv_t *pcvs_pcv;
    
    if (pcvs->nmax<=pcvs->n) {
        pcvs->nmax+=256;
        if (!(pcvs_pcv=(pcv_t *)realloc(pcvs->pcv,sizeof(pcv_t)*pcvs->nmax))) {
            trace(1,"addpcv: memory allocation error\n");
            free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
            return;
        }
        pcvs->pcv=pcvs_pcv;
    }
    pcvs->pcv[pcvs->n++]=*pcv;
}
/* read ngs antenna parameter file -------------------------------------------*/
static int od_rtk_rtkcmn_readngspcv(const char *file, pcvs_t *pcvs)
{
    FILE *fp;
    static const pcv_t pcv0={0};
    pcv_t pcv;
    double neu[3];
    int n=0;
    char buff[256];
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"ngs pcv file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (strlen(buff)>=62&&buff[61]=='|') continue;
        
        if (buff[0]!=' ') n=0; /* start line */
        if (++n==1) {
            pcv=pcv0;
            strncpy(pcv.type,buff,61); pcv.type[61]='\0';
        }
        else if (n==2) {
            if (od_rtk_rtkcmn_decodef(buff,3,neu)<3) continue;
            pcv.off[0][0]=neu[1];
            pcv.off[0][1]=neu[0];
            pcv.off[0][2]=neu[2];
        }
        else if (n==3) od_rtk_rtkcmn_decodef(buff,10,pcv.var[0]);
        else if (n==4) od_rtk_rtkcmn_decodef(buff,9,pcv.var[0]+10);
        else if (n==5) {
            if (od_rtk_rtkcmn_decodef(buff,3,neu)<3) continue;;
            pcv.off[1][0]=neu[1];
            pcv.off[1][1]=neu[0];
            pcv.off[1][2]=neu[2];
        }
        else if (n==6) od_rtk_rtkcmn_decodef(buff,10,pcv.var[1]);
        else if (n==7) {
            od_rtk_rtkcmn_decodef(buff,9,pcv.var[1]+10);
            od_rtk_rtkcmn_addpcv(&pcv,pcvs);
        }
    }
    fclose(fp);
    
    return 1;
}
/* read antex file ----------------------------------------------------------*/
static int od_rtk_rtkcmn_readantex(const char *file, pcvs_t *pcvs)
{
    FILE *fp;
    static const pcv_t pcv0={0};
    pcv_t pcv;
    double neu[3];
    int i,f,freq=0,state=0,freqs[]={1,2,5,0};
    char buff[256];
    
    trace(3,"readantex: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"antex pcv file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (strlen(buff)<60||strstr(buff+60,"COMMENT")) continue;
        
        if (strstr(buff+60,"START OF ANTENNA")) {
            pcv=pcv0;
            state=1;
        }
        if (strstr(buff+60,"END OF ANTENNA")) {
            od_rtk_rtkcmn_addpcv(&pcv,pcvs);
            state=0;
        }
        if (!state) continue;
        
        if (strstr(buff+60,"TYPE / SERIAL NO")) {
            strncpy(pcv.type,buff   ,20); pcv.type[20]='\0';
            strncpy(pcv.code,buff+20,20); pcv.code[20]='\0';
            if (!strncmp(pcv.code+3,"        ",8)) {
                pcv.sat=satid2no(pcv.code);
            }
        }
        else if (strstr(buff+60,"VALID FROM")) {
            if (!str2time(buff,0,43,&pcv.ts)) continue;
        }
        else if (strstr(buff+60,"VALID UNTIL")) {
            if (!str2time(buff,0,43,&pcv.te)) continue;
        }
        else if (strstr(buff+60,"START OF FREQUENCY")) {
            if (!pcv.sat&&buff[3]!='G') continue; /* only read rec ant for GPS */
            if (sscanf(buff+4,"%d",&f)<1) continue;
            for (i=0;freqs[i];i++) if (freqs[i]==f) break;
            if (freqs[i]) freq=i+1;
        }
        else if (strstr(buff+60,"END OF FREQUENCY")) {
            freq=0;
        }
        else if (strstr(buff+60,"NORTH / EAST / UP")) {
            if (freq<1||NFREQ<freq) continue;
            if (od_rtk_rtkcmn_decodef(buff,3,neu)<3) continue;
            pcv.off[freq-1][0]=neu[pcv.sat?0:1]; /* x or e */
            pcv.off[freq-1][1]=neu[pcv.sat?1:0]; /* y or n */
            pcv.off[freq-1][2]=neu[2];           /* z or u */
        }
        else if (strstr(buff,"NOAZI")) {
            if (freq<1||NFREQ<freq) continue;
            if ((i=od_rtk_rtkcmn_decodef(buff+8,19,pcv.var[freq-1]))<=0) continue;
            for (;i<19;i++) pcv.var[freq-1][i]=pcv.var[freq-1][i-1];
        }
    }
    fclose(fp);
    
    return 1;
}
/* read antenna parameters ------------------------------------------------------
* read antenna parameters
* args   : char   *file       I   antenna parameter file (antex)
*          pcvs_t *pcvs       IO  antenna parameters
* return : status (1:ok,0:file open error)
* notes  : file with the externsion .atx or .ATX is recognized as antex
*          file except for antex is recognized ngs antenna parameters
*          see reference [3]
*          only support non-azimuth-depedent parameters
*-----------------------------------------------------------------------------*/
extern int readpcv(const char *file, pcvs_t *pcvs)
{
    pcv_t *pcv;
    char *ext;
    int i,j,stat;
    
    trace(3,"readpcv: file=%s\n",file);
    
    if (!(ext=strrchr(file,'.'))) ext="";
    
    if (!strcmp(ext,".atx")||!strcmp(ext,".ATX")) {
        stat=od_rtk_rtkcmn_readantex(file,pcvs);
    }
    else {
        stat=od_rtk_rtkcmn_readngspcv(file,pcvs);
    }
    for (i=0;i<pcvs->n;i++) {
        pcv=pcvs->pcv+i;
        trace(4,"sat=%2d type=%20s code=%s off=%8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f\n",
              pcv->sat,pcv->type,pcv->code,pcv->off[0][0],pcv->off[0][1],
              pcv->off[0][2],pcv->off[1][0],pcv->off[1][1],pcv->off[1][2]);
        
        /* apply L2 to L3,L4,... if no pcv data */
        for (j=2;j<NFREQ;j++) { /* L3,L4,... */
            if (norm(pcv->off[j],3)>0.0) continue;
            matcpy(pcv->off[j],pcv->off[1], 3,1);
            matcpy(pcv->var[j],pcv->var[1],19,1);
        }
    }
    return stat;
}
/* search antenna parameter ----------------------------------------------------
* read satellite antenna phase center position
* args   : int    sat         I   satellite number (0: receiver antenna)
*          char   *type       I   antenna type for receiver antenna
*          gtime_t time       I   time to search parameters
*          pcvs_t *pcvs       IO  antenna parameters
* return : antenna parameter (NULL: no antenna)
*-----------------------------------------------------------------------------*/
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time,
                        const pcvs_t *pcvs)
{
    pcv_t *pcv;
    char buff[MAXANT],*types[2],*p;
    int i,j,n=0;
    
    trace(3,"searchpcv: sat=%2d type=%s\n",sat,type);
    
    if (sat) { /* search satellite antenna */
        for (i=0;i<pcvs->n;i++) {
            pcv=pcvs->pcv+i;
            if (pcv->sat!=sat) continue;
            if (pcv->ts.time!=0&&timediff(pcv->ts,time)>0.0) continue;
            if (pcv->te.time!=0&&timediff(pcv->te,time)<0.0) continue;
            return pcv;
        }
    }
    else {
        strcpy(buff,type);
        for (p=strtok(buff," ");p&&n<2;p=strtok(NULL," ")) types[n++]=p;
        if (n<=0) return NULL;
        
        /* search receiver antenna with radome at first */
        for (i=0;i<pcvs->n;i++) {
            pcv=pcvs->pcv+i;
            for (j=0;j<n;j++) if (!strstr(pcv->type,types[j])) break;
            if (j>=n) return pcv;
        }
        /* search receiver antenna without radome */
        for (i=0;i<pcvs->n;i++) {
            pcv=pcvs->pcv+i;
            if (strstr(pcv->type,types[0])!=pcv->type) continue;
            
            trace(2,"pcv without radome is used type=%s\n",type);
            return pcv;
        }
    }
    return NULL;
}
/* read station positions ------------------------------------------------------
* read positions from station position file
* args   : char  *file      I   station position file containing
*                               lat(deg) lon(deg) height(m) name in a line
*          char  *rcvs      I   station name
*          double *pos      O   station position {lat,lon,h} (rad/m)
*                               (all 0 if search error)
* return : none
*-----------------------------------------------------------------------------*/
extern void readpos(const char *file, const char *rcv, double *pos)
{
    static double poss[2048][3];
    static char stas[2048][16];
    FILE *fp;
    int i,j,len,np=0;
    char buff[256],str[256];
    
    trace(3,"readpos: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        fprintf(stderr,"reference position file open error : %s\n",file);
        return;
    }
    while (np<2048&&fgets(buff,sizeof(buff),fp)) {
        if (buff[0]=='%'||buff[0]=='#') continue;
        if (sscanf(buff,"%lf %lf %lf %s",&poss[np][0],&poss[np][1],&poss[np][2],
                   str)<4) continue;
        sprintf(stas[np++],"%.15s",str);
    }
    fclose(fp);
    len=(int)strlen(rcv);
    for (i=0;i<np;i++) {
        if (strncmp(stas[i],rcv,len)) continue;
        for (j=0;j<3;j++) pos[j]=poss[i][j];
        pos[0]*=D2R; pos[1]*=D2R;
        return;
    }
    pos[0]=pos[1]=pos[2]=0.0;
}
/* read blq record -----------------------------------------------------------*/
static int od_rtk_rtkcmn_readblqrecord(FILE *fp, double *odisp)
{
    double v[11];
    char buff[256];
    int i,n=0;
    
    while (fgets(buff,sizeof(buff),fp)) {
        if (!strncmp(buff,"$$",2)) continue;
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10)<11) continue;
        for (i=0;i<11;i++) odisp[n+i*6]=v[i];
        if (++n==6) return 1;
    }
    return 0;
}
/* read blq ocean tide loading parameters --------------------------------------
* read blq ocean tide loading parameters
* args   : char   *file       I   BLQ ocean tide loading parameter file
*          char   *sta        I   station name
*          double *odisp      O   ocean tide loading parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readblq(const char *file, const char *sta, double *odisp)
{
    FILE *fp;
    char buff[256],staname[32]="",name[32],*p;
    
    /* station name to upper case */
    sscanf(sta,"%16s",staname);
    for (p=staname;(*p=(char)toupper((int)(*p)));p++) ;
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"blq file open error: file=%s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (!strncmp(buff,"$$",2)||strlen(buff)<2) continue;
        
        if (sscanf(buff+2,"%16s",name)<1) continue;
        for (p=name;(*p=(char)toupper((int)(*p)));p++) ;
        if (strcmp(name,staname)) continue;
        
        /* read blq record */
        if (od_rtk_rtkcmn_readblqrecord(fp,odisp)) {
            fclose(fp);
            return 1;
        }
    }
    fclose(fp);
    trace(2,"no otl parameters: sta=%s file=%s\n",sta,file);
    return 0;
}
/* read earth rotation parameters ----------------------------------------------
* read earth rotation parameters
* args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readerp(const char *file, erp_t *erp)
{
    FILE *fp;
    erpd_t *erp_data;
    double v[14]={0};
    char buff[256];
    
    trace(3,"readerp: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"erp file open error: file=%s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                   v,v+1,v+2,v+3,v+4,v+5,v+6,v+7,v+8,v+9,v+10,v+11,v+12,v+13)<5) {
            continue;
        }
        if (erp->n>=erp->nmax) {
            erp->nmax=erp->nmax<=0?128:erp->nmax*2;
            erp_data=(erpd_t *)realloc(erp->data,sizeof(erpd_t)*erp->nmax);
            if (!erp_data) {
                free(erp->data); erp->data=NULL; erp->n=erp->nmax=0;
                fclose(fp);
                return 0;
            }
            erp->data=erp_data;
        }
        erp->data[erp->n].mjd=v[0];
        erp->data[erp->n].xp=v[1]*1E-6*AS2R;
        erp->data[erp->n].yp=v[2]*1E-6*AS2R;
        erp->data[erp->n].ut1_utc=v[3]*1E-7;
        erp->data[erp->n].lod=v[4]*1E-7;
        erp->data[erp->n].xpr=v[12]*1E-6*AS2R;
        erp->data[erp->n++].ypr=v[13]*1E-6*AS2R;
    }
    fclose(fp);
    return 1;
}
/* get earth rotation parameter values -----------------------------------------
* get earth rotation parameter values
* args   : erp_t  *erp        I   earth rotation parameters
*          gtime_t time       I   time (gpst)
*          double *erpv       O   erp values {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int geterp(const erp_t *erp, gtime_t time, double *erpv)
{
    const double ep[]={2000,1,1,12,0,0};
    double mjd,day,a;
    int i,j,k;
    
    trace(4,"geterp:\n");
    
    if (erp->n<=0) return 0;
    
    mjd=51544.5+(timediff(gpst2utc(time),epoch2time(ep)))/86400.0;
    
    if (mjd<=erp->data[0].mjd) {
        day=mjd-erp->data[0].mjd;
        erpv[0]=erp->data[0].xp     +erp->data[0].xpr*day;
        erpv[1]=erp->data[0].yp     +erp->data[0].ypr*day;
        erpv[2]=erp->data[0].ut1_utc-erp->data[0].lod*day;
        erpv[3]=erp->data[0].lod;
        return 1;
    }
    if (mjd>=erp->data[erp->n-1].mjd) {
        day=mjd-erp->data[erp->n-1].mjd;
        erpv[0]=erp->data[erp->n-1].xp     +erp->data[erp->n-1].xpr*day;
        erpv[1]=erp->data[erp->n-1].yp     +erp->data[erp->n-1].ypr*day;
        erpv[2]=erp->data[erp->n-1].ut1_utc-erp->data[erp->n-1].lod*day;
        erpv[3]=erp->data[erp->n-1].lod;
        return 1;
    }
    for (j=0,k=erp->n-1;j<k-1;) {
        i=(j+k)/2;
        if (mjd<erp->data[i].mjd) k=i; else j=i;
    }
    if (erp->data[j].mjd==erp->data[j+1].mjd) {
        a=0.5;
    }
    else {
        a=(mjd-erp->data[j].mjd)/(erp->data[j+1].mjd-erp->data[j].mjd);
    }
    erpv[0]=(1.0-a)*erp->data[j].xp     +a*erp->data[j+1].xp;
    erpv[1]=(1.0-a)*erp->data[j].yp     +a*erp->data[j+1].yp;
    erpv[2]=(1.0-a)*erp->data[j].ut1_utc+a*erp->data[j+1].ut1_utc;
    erpv[3]=(1.0-a)*erp->data[j].lod    +a*erp->data[j+1].lod;
    return 1;
}
/* compare ephemeris ---------------------------------------------------------*/
static int od_rtk_rtkcmn_cmpeph(const void *p1, const void *p2)
{
    eph_t *q1=(eph_t *)p1,*q2=(eph_t *)p2;
    return q1->ttr.time!=q2->ttr.time?(int)(q1->ttr.time-q2->ttr.time):
           (q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
            q1->sat-q2->sat);
}
/* sort and unique ephemeris -------------------------------------------------*/
static void od_rtk_rtkcmn_uniqeph(nav_t *nav)
{
    eph_t *nav_eph;
    int i,j;
    
    trace(3,"uniqeph: n=%d\n",nav->n);
    
    if (nav->n<=0) return;
    
    qsort(nav->eph,nav->n,sizeof(eph_t),od_rtk_rtkcmn_cmpeph);
    
    for (i=1,j=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=nav->eph[j].sat||
            nav->eph[i].iode!=nav->eph[j].iode) {
            nav->eph[++j]=nav->eph[i];
        }
    }
    nav->n=j+1;
    
    if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->n))) {
        trace(1,"uniqeph malloc error n=%d\n",nav->n);
        free(nav->eph); nav->eph=NULL; nav->n=nav->nmax=0;
        return;
    }
    nav->eph=nav_eph;
    nav->nmax=nav->n;
    
    trace(4,"uniqeph: n=%d\n",nav->n);
}
/* compare glonass ephemeris -------------------------------------------------*/
static int od_rtk_rtkcmn_cmpgeph(const void *p1, const void *p2)
{
    geph_t *q1=(geph_t *)p1,*q2=(geph_t *)p2;
    return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
           (q1->toe.time!=q2->toe.time?(int)(q1->toe.time-q2->toe.time):
            q1->sat-q2->sat);
}
/* sort and unique glonass ephemeris -----------------------------------------*/
static void od_rtk_rtkcmn_uniqgeph(nav_t *nav)
{
    geph_t *nav_geph;
    int i,j;
    
    trace(3,"uniqgeph: ng=%d\n",nav->ng);
    
    if (nav->ng<=0) return;
    
    qsort(nav->geph,nav->ng,sizeof(geph_t),od_rtk_rtkcmn_cmpgeph);
    
    for (i=j=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=nav->geph[j].sat||
            nav->geph[i].toe.time!=nav->geph[j].toe.time||
            nav->geph[i].svh!=nav->geph[j].svh) {
            nav->geph[++j]=nav->geph[i];
        }
    }
    nav->ng=j+1;
    
    if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ng))) {
        trace(1,"uniqgeph malloc error ng=%d\n",nav->ng);
        free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
        return;
    }
    nav->geph=nav_geph;
    nav->ngmax=nav->ng;
    
    trace(4,"uniqgeph: ng=%d\n",nav->ng);
}
/* compare sbas ephemeris ----------------------------------------------------*/
static int od_rtk_rtkcmn_cmpseph(const void *p1, const void *p2)
{
    seph_t *q1=(seph_t *)p1,*q2=(seph_t *)p2;
    return q1->tof.time!=q2->tof.time?(int)(q1->tof.time-q2->tof.time):
           (q1->t0.time!=q2->t0.time?(int)(q1->t0.time-q2->t0.time):
            q1->sat-q2->sat);
}
/* sort and unique sbas ephemeris --------------------------------------------*/
static void od_rtk_rtkcmn_uniqseph(nav_t *nav)
{
    seph_t *nav_seph;
    int i,j;
    
    trace(3,"uniqseph: ns=%d\n",nav->ns);
    
    if (nav->ns<=0) return;
    
    qsort(nav->seph,nav->ns,sizeof(seph_t),od_rtk_rtkcmn_cmpseph);
    
    for (i=j=0;i<nav->ns;i++) {
        if (nav->seph[i].sat!=nav->seph[j].sat||
            nav->seph[i].t0.time!=nav->seph[j].t0.time) {
            nav->seph[++j]=nav->seph[i];
        }
    }
    nav->ns=j+1;
    
    if (!(nav_seph=(seph_t *)realloc(nav->seph,sizeof(seph_t)*nav->ns))) {
        trace(1,"uniqseph malloc error ns=%d\n",nav->ns);
        free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
        return;
    }
    nav->seph=nav_seph;
    nav->nsmax=nav->ns;
    
    trace(4,"uniqseph: ns=%d\n",nav->ns);
}
/* unique ephemerides ----------------------------------------------------------
* unique ephemerides in navigation data
* args   : nav_t *nav    IO     navigation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern void uniqnav(nav_t *nav)
{
    trace(3,"uniqnav: neph=%d ngeph=%d nseph=%d\n",nav->n,nav->ng,nav->ns);
    
    /* unique ephemeris */
    od_rtk_rtkcmn_uniqeph (nav);
    od_rtk_rtkcmn_uniqgeph(nav);
    od_rtk_rtkcmn_uniqseph(nav);
}
/* compare observation data -------------------------------------------------*/
static int od_rtk_rtkcmn_cmpobs(const void *p1, const void *p2)
{
    obsd_t *q1=(obsd_t *)p1,*q2=(obsd_t *)p2;
    double tt=timediff(q1->time,q2->time);
    if (fabs(tt)>DTTOL) return tt<0?-1:1;
    if (q1->rcv!=q2->rcv) return (int)q1->rcv-(int)q2->rcv;
    return (int)q1->sat-(int)q2->sat;
}
/* sort and unique observation data --------------------------------------------
* sort and unique observation data by time, rcv, sat
* args   : obs_t *obs    IO     observation data
* return : number of epochs
*-----------------------------------------------------------------------------*/
extern int sortobs(obs_t *obs)
{
    int i,j,n;
    
    trace(3,"sortobs: nobs=%d\n",obs->n);
    
    if (obs->n<=0) return 0;
    
    qsort(obs->data,obs->n,sizeof(obsd_t),od_rtk_rtkcmn_cmpobs);
    
    /* delete duplicated data */
    for (i=j=0;i<obs->n;i++) {
        if (obs->data[i].sat!=obs->data[j].sat||
            obs->data[i].rcv!=obs->data[j].rcv||
            timediff(obs->data[i].time,obs->data[j].time)!=0.0) {
            obs->data[++j]=obs->data[i];
        }
    }
    obs->n=j+1;
    
    for (i=n=0;i<obs->n;i=j,n++) {
        for (j=i+1;j<obs->n;j++) {
            if (timediff(obs->data[j].time,obs->data[i].time)>DTTOL) break;
        }
    }
    return n;
}
/* screen by time --------------------------------------------------------------
* screening by time start, time end, and time interval
* args   : gtime_t time  I      time
*          gtime_t ts    I      time start (ts.time==0:no screening by ts)
*          gtime_t te    I      time end   (te.time==0:no screening by te)
*          double  tint  I      time interval (s) (0.0:no screen by tint)
* return : 1:on condition, 0:not on condition
*-----------------------------------------------------------------------------*/
extern int screent(gtime_t time, gtime_t ts, gtime_t te, double tint)
{
    return (tint<=0.0||fmod(time2gpst(time,NULL)+DTTOL,tint)<=DTTOL*2.0)&&
           (ts.time==0||timediff(time,ts)>=-DTTOL)&&
           (te.time==0||timediff(time,te)<  DTTOL);
}
/* read/save navigation data ---------------------------------------------------
* save or load navigation data
* args   : char    file  I      file path
*          nav_t   nav   O/I    navigation data
* return : status (1:ok,0:no file)
*-----------------------------------------------------------------------------*/
extern int readnav(const char *file, nav_t *nav)
{
    FILE *fp;
    eph_t eph0={0};
    geph_t geph0={0};
    char buff[4096],*p;
    long toe_time,tof_time,toc_time,ttr_time;
    int i,sat,prn;
    
    trace(3,"loadnav: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) return 0;
    
    while (fgets(buff,sizeof(buff),fp)) {
        if (!strncmp(buff,"IONUTC",6)) {
            for (i=0;i<8;i++) nav->ion_gps[i]=0.0;
            for (i=0;i<8;i++) nav->utc_gps[i]=0.0;
            sscanf(buff,"IONUTC,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &nav->ion_gps[0],&nav->ion_gps[1],&nav->ion_gps[2],&nav->ion_gps[3],
                   &nav->ion_gps[4],&nav->ion_gps[5],&nav->ion_gps[6],&nav->ion_gps[7],
                   &nav->utc_gps[0],&nav->utc_gps[1],&nav->utc_gps[2],&nav->utc_gps[3],
                   &nav->utc_gps[4]);
            continue;   
        }
        if ((p=strchr(buff,','))) *p='\0'; else continue;
        if (!(sat=satid2no(buff))) continue;
        if (satsys(sat,&prn)==SYS_GLO) {
            nav->geph[prn-1]=geph0;
            nav->geph[prn-1].sat=sat;
            toe_time=tof_time=0;
            sscanf(p+1,"%d,%d,%d,%d,%d,%ld,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,"
                        "%lf,%lf,%lf,%lf",
                   &nav->geph[prn-1].iode,&nav->geph[prn-1].frq,&nav->geph[prn-1].svh,
                   &nav->geph[prn-1].sva,&nav->geph[prn-1].age,
                   &toe_time,&tof_time,
                   &nav->geph[prn-1].pos[0],&nav->geph[prn-1].pos[1],&nav->geph[prn-1].pos[2],
                   &nav->geph[prn-1].vel[0],&nav->geph[prn-1].vel[1],&nav->geph[prn-1].vel[2],
                   &nav->geph[prn-1].acc[0],&nav->geph[prn-1].acc[1],&nav->geph[prn-1].acc[2],
                   &nav->geph[prn-1].taun  ,&nav->geph[prn-1].gamn  ,&nav->geph[prn-1].dtaun);
            nav->geph[prn-1].toe.time=toe_time;
            nav->geph[prn-1].tof.time=tof_time;
        }
        else {
            nav->eph[sat-1]=eph0;
            nav->eph[sat-1].sat=sat;
            toe_time=toc_time=ttr_time=0;
            sscanf(p+1,"%d,%d,%d,%d,%ld,%ld,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,"
                        "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d",
                   &nav->eph[sat-1].iode,&nav->eph[sat-1].iodc,&nav->eph[sat-1].sva ,
                   &nav->eph[sat-1].svh ,
                   &toe_time,&toc_time,&ttr_time,
                   &nav->eph[sat-1].A   ,&nav->eph[sat-1].e   ,&nav->eph[sat-1].i0  ,
                   &nav->eph[sat-1].OMG0,&nav->eph[sat-1].omg ,&nav->eph[sat-1].M0  ,
                   &nav->eph[sat-1].deln,&nav->eph[sat-1].OMGd,&nav->eph[sat-1].idot,
                   &nav->eph[sat-1].crc ,&nav->eph[sat-1].crs ,&nav->eph[sat-1].cuc ,
                   &nav->eph[sat-1].cus ,&nav->eph[sat-1].cic ,&nav->eph[sat-1].cis ,
                   &nav->eph[sat-1].toes,&nav->eph[sat-1].fit ,&nav->eph[sat-1].f0  ,
                   &nav->eph[sat-1].f1  ,&nav->eph[sat-1].f2  ,&nav->eph[sat-1].tgd[0],
                   &nav->eph[sat-1].code, &nav->eph[sat-1].flag);
            nav->eph[sat-1].toe.time=toe_time;
            nav->eph[sat-1].toc.time=toc_time;
            nav->eph[sat-1].ttr.time=ttr_time;
        }
    }
    fclose(fp);
    return 1;
}
extern int savenav(const char *file, const nav_t *nav)
{
    FILE *fp;
    int i;
    char id[32];
    
    trace(3,"savenav: file=%s\n",file);
    
    if (!(fp=fopen(file,"w"))) return 0;
    
    for (i=0;i<MAXSAT;i++) {
        if (nav->eph[i].ttr.time==0) continue;
        satno2id(nav->eph[i].sat,id);
        fprintf(fp,"%s,%d,%d,%d,%d,%d,%d,%d,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
                   "%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
                   "%.14E,%.14E,%.14E,%.14E,%.14E,%d,%d\n",
                id,nav->eph[i].iode,nav->eph[i].iodc,nav->eph[i].sva ,
                nav->eph[i].svh ,(int)nav->eph[i].toe.time,
                (int)nav->eph[i].toc.time,(int)nav->eph[i].ttr.time,
                nav->eph[i].A   ,nav->eph[i].e  ,nav->eph[i].i0  ,nav->eph[i].OMG0,
                nav->eph[i].omg ,nav->eph[i].M0 ,nav->eph[i].deln,nav->eph[i].OMGd,
                nav->eph[i].idot,nav->eph[i].crc,nav->eph[i].crs ,nav->eph[i].cuc ,
                nav->eph[i].cus ,nav->eph[i].cic,nav->eph[i].cis ,nav->eph[i].toes,
                nav->eph[i].fit ,nav->eph[i].f0 ,nav->eph[i].f1  ,nav->eph[i].f2  ,
                nav->eph[i].tgd[0],nav->eph[i].code,nav->eph[i].flag);
    }
    for (i=0;i<MAXPRNGLO;i++) {
        if (nav->geph[i].tof.time==0) continue;
        satno2id(nav->geph[i].sat,id);
        fprintf(fp,"%s,%d,%d,%d,%d,%d,%d,%d,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
                   "%.14E,%.14E,%.14E,%.14E,%.14E,%.14E\n",
                id,nav->geph[i].iode,nav->geph[i].frq,nav->geph[i].svh,
                nav->geph[i].sva,nav->geph[i].age,(int)nav->geph[i].toe.time,
                (int)nav->geph[i].tof.time,
                nav->geph[i].pos[0],nav->geph[i].pos[1],nav->geph[i].pos[2],
                nav->geph[i].vel[0],nav->geph[i].vel[1],nav->geph[i].vel[2],
                nav->geph[i].acc[0],nav->geph[i].acc[1],nav->geph[i].acc[2],
                nav->geph[i].taun,nav->geph[i].gamn,nav->geph[i].dtaun);
    }
    fprintf(fp,"IONUTC,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,%.14E,"
               "%.14E,%.14E,%.14E,%.0f",
            nav->ion_gps[0],nav->ion_gps[1],nav->ion_gps[2],nav->ion_gps[3],
            nav->ion_gps[4],nav->ion_gps[5],nav->ion_gps[6],nav->ion_gps[7],
            nav->utc_gps[0],nav->utc_gps[1],nav->utc_gps[2],nav->utc_gps[3],
            nav->utc_gps[4]);
    
    fclose(fp);
    return 1;
}
/* free observation data -------------------------------------------------------
* free memory for observation data
* args   : obs_t *obs    IO     observation data
* return : none
*-----------------------------------------------------------------------------*/
extern void freeobs(obs_t *obs)
{
    free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
}
/* free navigation data ---------------------------------------------------------
* free memory for navigation data
* args   : nav_t *nav    IO     navigation data
*          int   opt     I      option (or of followings)
*                               (0x01: gps/qzs ephmeris, 0x02: glonass ephemeris,
*                                0x04: sbas ephemeris,   0x08: precise ephemeris,
*                                0x10: precise clock     0x20: almanac,
*                                0x40: tec data)
* return : none
*-----------------------------------------------------------------------------*/
extern void freenav(nav_t *nav, int opt)
{
    if (opt&0x01) {free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;}
    if (opt&0x02) {free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;}
    if (opt&0x04) {free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;}
    if (opt&0x08) {free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;}
    if (opt&0x10) {free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;}
    if (opt&0x20) {free(nav->alm ); nav->alm =NULL; nav->na=nav->namax=0;}
    if (opt&0x40) {free(nav->tec ); nav->tec =NULL; nav->nt=nav->ntmax=0;}
}
/* debug trace functions -----------------------------------------------------*/
#ifdef TRACE

static FILE *od_rtk_rtkcmn_fp_trace=NULL;     /* file pointer of trace */
static char od_rtk_rtkcmn_file_trace[1024];   /* trace file */
static int od_rtk_rtkcmn_level_trace=0;       /* level of trace */
static uint32_t od_rtk_rtkcmn_tick_trace=0;   /* tick time at traceopen (ms) */
static gtime_t od_rtk_rtkcmn_time_trace={0};  /* time at traceopen */
static lock_t od_rtk_rtkcmn_lock_trace;       /* lock for trace */

static void od_rtk_rtkcmn_traceswap(void)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];
    
    lock(&od_rtk_rtkcmn_lock_trace);
    
    if ((int)(time2gpst(time      ,NULL)/INT_SWAP_TRAC)==
        (int)(time2gpst(od_rtk_rtkcmn_time_trace,NULL)/INT_SWAP_TRAC)) {
        unlock(&od_rtk_rtkcmn_lock_trace);
        return;
    }
    od_rtk_rtkcmn_time_trace=time;
    
    if (!reppath(od_rtk_rtkcmn_file_trace,path,time,"","")) {
        unlock(&od_rtk_rtkcmn_lock_trace);
        return;
    }
    if (od_rtk_rtkcmn_fp_trace) fclose(od_rtk_rtkcmn_fp_trace);
    
    if (!(od_rtk_rtkcmn_fp_trace=fopen(path,"w"))) {
        od_rtk_rtkcmn_fp_trace=stderr;
    }
    unlock(&od_rtk_rtkcmn_lock_trace);
}
extern void traceopen(const char *file)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];
    
    reppath(file,path,time,"","");
    if (!*path||!(od_rtk_rtkcmn_fp_trace=fopen(path,"w"))) od_rtk_rtkcmn_fp_trace=stderr;
    strcpy(od_rtk_rtkcmn_file_trace,file);
    od_rtk_rtkcmn_tick_trace=tickget();
    od_rtk_rtkcmn_time_trace=time;
    initlock(&od_rtk_rtkcmn_lock_trace);
}
extern void traceclose(void)
{
    if (od_rtk_rtkcmn_fp_trace&&od_rtk_rtkcmn_fp_trace!=stderr) fclose(od_rtk_rtkcmn_fp_trace);
    od_rtk_rtkcmn_fp_trace=NULL;
    od_rtk_rtkcmn_file_trace[0]='\0';
}
extern void tracelevel(int level)
{
    od_rtk_rtkcmn_level_trace=level;
}
extern void trace(int level, const char *format, ...)
{
    va_list ap;
    
    /* print error message to stderr */
    if (level<=1) {
        va_start(ap,format); vfprintf(stderr,format,ap); va_end(ap);
    }
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    od_rtk_rtkcmn_traceswap();
    fprintf(od_rtk_rtkcmn_fp_trace,"%d ",level);
    va_start(ap,format); vfprintf(od_rtk_rtkcmn_fp_trace,format,ap); va_end(ap);
    fflush(od_rtk_rtkcmn_fp_trace);
}
extern void tracet(int level, const char *format, ...)
{
    va_list ap;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    od_rtk_rtkcmn_traceswap();
    fprintf(od_rtk_rtkcmn_fp_trace,"%d %9.3f: ",level,(tickget()-od_rtk_rtkcmn_tick_trace)/1000.0);
    va_start(ap,format); vfprintf(od_rtk_rtkcmn_fp_trace,format,ap); va_end(ap);
    fflush(od_rtk_rtkcmn_fp_trace);
}
extern void tracemat(int level, const double *A, int n, int m, int p, int q)
{
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    matfprint(A,n,m,p,q,od_rtk_rtkcmn_fp_trace); fflush(od_rtk_rtkcmn_fp_trace);
}
extern void traceobs(int level, const obsd_t *obs, int n)
{
    char str[64],id[16];
    int i;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    for (i=0;i<n;i++) {
        time2str(obs[i].time,str,3);
        satno2id(obs[i].sat,id);
        fprintf(od_rtk_rtkcmn_fp_trace," (%2d) %s %-3s rcv%d %13.3f %13.3f %13.3f %13.3f %d %d %d %d %3.1f %3.1f\n",
              i+1,str,id,obs[i].rcv,obs[i].L[0],obs[i].L[1],obs[i].P[0],
              obs[i].P[1],obs[i].LLI[0],obs[i].LLI[1],obs[i].code[0],
              obs[i].code[1],obs[i].SNR[0]*SNR_UNIT,obs[i].SNR[1]*SNR_UNIT);
    }
    fflush(od_rtk_rtkcmn_fp_trace);
}
extern void tracenav(int level, const nav_t *nav)
{
    char s1[64],s2[64],id[16];
    int i;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    for (i=0;i<nav->n;i++) {
        time2str(nav->eph[i].toe,s1,0);
        time2str(nav->eph[i].ttr,s2,0);
        satno2id(nav->eph[i].sat,id);
        fprintf(od_rtk_rtkcmn_fp_trace,"(%3d) %-3s : %s %s %3d %3d %02x\n",i+1,
                id,s1,s2,nav->eph[i].iode,nav->eph[i].iodc,nav->eph[i].svh);
    }
    fprintf(od_rtk_rtkcmn_fp_trace,"(ion) %9.4e %9.4e %9.4e %9.4e\n",nav->ion_gps[0],
            nav->ion_gps[1],nav->ion_gps[2],nav->ion_gps[3]);
    fprintf(od_rtk_rtkcmn_fp_trace,"(ion) %9.4e %9.4e %9.4e %9.4e\n",nav->ion_gps[4],
            nav->ion_gps[5],nav->ion_gps[6],nav->ion_gps[7]);
    fprintf(od_rtk_rtkcmn_fp_trace,"(ion) %9.4e %9.4e %9.4e %9.4e\n",nav->ion_gal[0],
            nav->ion_gal[1],nav->ion_gal[2],nav->ion_gal[3]);
}
extern void tracegnav(int level, const nav_t *nav)
{
    char s1[64],s2[64],id[16];
    int i;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    for (i=0;i<nav->ng;i++) {
        time2str(nav->geph[i].toe,s1,0);
        time2str(nav->geph[i].tof,s2,0);
        satno2id(nav->geph[i].sat,id);
        fprintf(od_rtk_rtkcmn_fp_trace,"(%3d) %-3s : %s %s %2d %2d %8.3f\n",i+1,
                id,s1,s2,nav->geph[i].frq,nav->geph[i].svh,nav->geph[i].taun*1E6);
    }
}
extern void tracehnav(int level, const nav_t *nav)
{
    char s1[64],s2[64],id[16];
    int i;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    for (i=0;i<nav->ns;i++) {
        time2str(nav->seph[i].t0,s1,0);
        time2str(nav->seph[i].tof,s2,0);
        satno2id(nav->seph[i].sat,id);
        fprintf(od_rtk_rtkcmn_fp_trace,"(%3d) %-3s : %s %s %2d %2d\n",i+1,
                id,s1,s2,nav->seph[i].svh,nav->seph[i].sva);
    }
}
extern void tracepeph(int level, const nav_t *nav)
{
    char s[64],id[16];
    int i,j;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    
    for (i=0;i<nav->ne;i++) {
        time2str(nav->peph[i].time,s,0);
        for (j=0;j<MAXSAT;j++) {
            satno2id(j+1,id);
            fprintf(od_rtk_rtkcmn_fp_trace,"%-3s %d %-3s %13.3f %13.3f %13.3f %13.3f %6.3f %6.3f %6.3f %6.3f\n",
                    s,nav->peph[i].index,id,
                    nav->peph[i].pos[j][0],nav->peph[i].pos[j][1],
                    nav->peph[i].pos[j][2],nav->peph[i].pos[j][3]*1E9,
                    nav->peph[i].std[j][0],nav->peph[i].std[j][1],
                    nav->peph[i].std[j][2],nav->peph[i].std[j][3]*1E9);
        }
    }
}
extern void tracepclk(int level, const nav_t *nav)
{
    char s[64],id[16];
    int i,j;
    
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    
    for (i=0;i<nav->nc;i++) {
        time2str(nav->pclk[i].time,s,0);
        for (j=0;j<MAXSAT;j++) {
            satno2id(j+1,id);
            fprintf(od_rtk_rtkcmn_fp_trace,"%-3s %d %-3s %13.3f %6.3f\n",
                    s,nav->pclk[i].index,id,
                    nav->pclk[i].clk[j][0]*1E9,nav->pclk[i].std[j][0]*1E9);
        }
    }
}
extern void traceb(int level, const uint8_t *p, int n)
{
    int i;
    if (!od_rtk_rtkcmn_fp_trace||level>od_rtk_rtkcmn_level_trace) return;
    for (i=0;i<n;i++) fprintf(od_rtk_rtkcmn_fp_trace,"%02X%s",*p++,i%8==7?" ":"");
    fprintf(od_rtk_rtkcmn_fp_trace,"\n");
}
#else
extern void traceopen(const char *file) {}
extern void traceclose(void) {}
extern void tracelevel(int level) {}
extern void trace   (int level, const char *format, ...) {}
extern void tracet  (int level, const char *format, ...) {}
extern void tracemat(int level, const double *A, int n, int m, int p, int q) {}
extern void traceobs(int level, const obsd_t *obs, int n) {}
extern void tracenav(int level, const nav_t *nav) {}
extern void tracegnav(int level, const nav_t *nav) {}
extern void tracehnav(int level, const nav_t *nav) {}
extern void tracepeph(int level, const nav_t *nav) {}
extern void tracepclk(int level, const nav_t *nav) {}
extern void traceb  (int level, const uint8_t *p, int n) {}

#endif /* TRACE */

/* execute command -------------------------------------------------------------
* execute command line by operating system shell
* args   : char   *cmd      I   command line
* return : execution status (0:ok,0>:error)
*-----------------------------------------------------------------------------*/
extern int execcmd(const char *cmd)
{
#ifdef WIN32
    PROCESS_INFORMATION info;
    STARTUPINFO si={0};
    DWORD stat;
    char cmds[1024];
    
    trace(3,"execcmd: cmd=%s\n",cmd);
    
    si.cb=sizeof(si);
    sprintf(cmds,"cmd /c %s",cmd);
    if (!CreateProcess(NULL,(LPTSTR)cmds,NULL,NULL,FALSE,CREATE_NO_WINDOW,NULL,
                       NULL,&si,&info)) return -1;
    WaitForSingleObject(info.hProcess,INFINITE);
    if (!GetExitCodeProcess(info.hProcess,&stat)) stat=-1;
    CloseHandle(info.hProcess);
    CloseHandle(info.hThread);
    return (int)stat;
#else
    trace(3,"execcmd: cmd=%s\n",cmd);
    
    return system(cmd);
#endif
}
/* expand file path ------------------------------------------------------------
* expand file path with wild-card (*) in file
* args   : char   *path     I   file path to expand (captal insensitive)
*          char   *paths    O   expanded file paths
*          int    nmax      I   max number of expanded file paths
* return : number of expanded file paths
* notes  : the order of expanded files is alphabetical order
*-----------------------------------------------------------------------------*/
extern int expath(const char *path, char *paths[], int nmax)
{
    int i,j,n=0;
    char tmp[1024];
#ifdef WIN32
    WIN32_FIND_DATA file;
    HANDLE h;
    char dir[1024]="",*p;
    
    trace(3,"expath  : path=%s nmax=%d\n",path,nmax);
    
    if ((p=strrchr(path,'\\'))) {
        strncpy(dir,path,p-path+1); dir[p-path+1]='\0';
    }
    if ((h=FindFirstFile((LPCTSTR)path,&file))==INVALID_HANDLE_VALUE) {
        strcpy(paths[0],path);
        return 1;
    }
    sprintf(paths[n++],"%s%s",dir,file.cFileName);
    while (FindNextFile(h,&file)&&n<nmax) {
        if (file.dwFileAttributes&FILE_ATTRIBUTE_DIRECTORY) continue;
        sprintf(paths[n++],"%s%s",dir,file.cFileName);
    }
    FindClose(h);
#else
    struct dirent *d;
    DIR *dp;
    const char *file=path;
    char dir[1024]="",s1[1024],s2[1024],*p,*q,*r;
    
    trace(3,"expath  : path=%s nmax=%d\n",path,nmax);
    
    if ((p=strrchr(path,'/'))||(p=strrchr(path,'\\'))) {
        file=p+1; strncpy(dir,path,p-path+1); dir[p-path+1]='\0';
    }
    if (!(dp=opendir(*dir?dir:"."))) return 0;
    while ((d=readdir(dp))) {
        if (*(d->d_name)=='.') continue;
        sprintf(s1,"^%s$",d->d_name);
        sprintf(s2,"^%s$",file);
        for (p=s1;*p;p++) *p=(char)tolower((int)*p);
        for (p=s2;*p;p++) *p=(char)tolower((int)*p);
        
        for (p=s1,q=strtok_r(s2,"*",&r);q;q=strtok_r(NULL,"*",&r)) {
            if ((p=strstr(p,q))) p+=strlen(q); else break;
        }
        if (p&&n<nmax) sprintf(paths[n++],"%s%s",dir,d->d_name);
    }
    closedir(dp);
#endif
    /* sort paths in alphabetical order */
    for (i=0;i<n-1;i++) {
        for (j=i+1;j<n;j++) {
            if (strcmp(paths[i],paths[j])>0) {
                strcpy(tmp,paths[i]);
                strcpy(paths[i],paths[j]);
                strcpy(paths[j],tmp);
            }
        }
    }
    for (i=0;i<n;i++) trace(3,"expath  : file=%s\n",paths[i]);
    
    return n;
}
/* generate local directory recursively --------------------------------------*/
static int od_rtk_rtkcmn_mkdir_r(const char *dir)
{
    char pdir[1024],*p;

#ifdef WIN32
    HANDLE h;
    WIN32_FIND_DATA data;
    
    if (!*dir||!strcmp(dir+1,":\\")) return 1;
    
    sprintf(pdir,"%.1023s",dir);
    if ((p=strrchr(pdir,FILEPATHSEP))) {
        *p='\0';
        h=FindFirstFile(pdir,&data);
        if (h==INVALID_HANDLE_VALUE) {
            if (!od_rtk_rtkcmn_mkdir_r(pdir)) return 0;
        }
        else FindClose(h);
    }
    if (CreateDirectory(dir,NULL)||GetLastError()==ERROR_ALREADY_EXISTS) {
        return 1;
    }
#else
    FILE *fp;
    
    if (!*dir) return 1;
    
    sprintf(pdir,"%.1023s",dir);
    if ((p=strrchr(pdir,FILEPATHSEP))) {
        *p='\0';
        if (!(fp=fopen(pdir,"r"))) {
            if (!od_rtk_rtkcmn_mkdir_r(pdir)) return 0;
        }
        else fclose(fp);
    }
    if (!mkdir(dir,0777)||errno==EEXIST) return 1;
#endif
    trace(2,"directory generation error: dir=%s\n",dir);
    return 0;
}
/* create directory ------------------------------------------------------------
* create directory if not exists
* args   : char   *path     I   file path to be saved
* return : none
* notes  : recursively.
*-----------------------------------------------------------------------------*/
extern void createdir(const char *path)
{
    char buff[1024],*p;
    
    tracet(3,"createdir: path=%s\n",path);
    
    strcpy(buff,path);
    if (!(p=strrchr(buff,FILEPATHSEP))) return;
    *p='\0';
    
    od_rtk_rtkcmn_mkdir_r(buff);
}
/* replace string ------------------------------------------------------------*/
static int od_rtk_rtkcmn_repstr(char *str, const char *pat, const char *rep)
{
    int len=(int)strlen(pat);
    char buff[1024],*p,*q,*r;
    
    for (p=str,r=buff;*p;p=q+len) {
        if (!(q=strstr(p,pat))) break;
        strncpy(r,p,q-p);
        r+=q-p;
        r+=sprintf(r,"%s",rep);
    }
    if (p<=str) return 0;
    strcpy(r,p);
    strcpy(str,buff);
    return 1;
}
/* replace keywords in file path -----------------------------------------------
* replace keywords in file path with date, time, rover and base station id
* args   : char   *path     I   file path (see below)
*          char   *rpath    O   file path in which keywords replaced (see below)
*          gtime_t time     I   time (gpst)  (time.time==0: not replaced)
*          char   *rov      I   rover id string        ("": not replaced)
*          char   *base     I   base station id string ("": not replaced)
* return : status (1:keywords replaced, 0:no valid keyword in the path,
*                  -1:no valid time)
* notes  : the following keywords in path are replaced by date, time and name
*              %Y -> yyyy : year (4 digits) (1900-2099)
*              %y -> yy   : year (2 digits) (00-99)
*              %m -> mm   : month           (01-12)
*              %d -> dd   : day of month    (01-31)
*              %h -> hh   : hours           (00-23)
*              %M -> mm   : minutes         (00-59)
*              %S -> ss   : seconds         (00-59)
*              %n -> ddd  : day of year     (001-366)
*              %W -> wwww : gps week        (0001-9999)
*              %D -> d    : day of gps week (0-6)
*              %H -> h    : hour code       (a=0,b=1,c=2,...,x=23)
*              %ha-> hh   : 3 hours         (00,03,06,...,21)
*              %hb-> hh   : 6 hours         (00,06,12,18)
*              %hc-> hh   : 12 hours        (00,12)
*              %t -> mm   : 15 minutes      (00,15,30,45)
*              %r -> rrrr : rover id
*              %b -> bbbb : base station id
*-----------------------------------------------------------------------------*/
extern int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
                   const char *base)
{
    double ep[6],ep0[6]={2000,1,1,0,0,0};
    int week,dow,doy,stat=0;
    char rep[64];
    
    strcpy(rpath,path);
    
    if (!strstr(rpath,"%")) return 0;
    if (*rov ) stat|=od_rtk_rtkcmn_repstr(rpath,"%r",rov );
    if (*base) stat|=od_rtk_rtkcmn_repstr(rpath,"%b",base);
    if (time.time!=0) {
        time2epoch(time,ep);
        ep0[0]=ep[0];
        dow=(int)floor(time2gpst(time,&week)/86400.0);
        doy=(int)floor(timediff(time,epoch2time(ep0))/86400.0)+1;
        sprintf(rep,"%02d",  ((int)ep[3]/3)*3);   stat|=od_rtk_rtkcmn_repstr(rpath,"%ha",rep);
        sprintf(rep,"%02d",  ((int)ep[3]/6)*6);   stat|=od_rtk_rtkcmn_repstr(rpath,"%hb",rep);
        sprintf(rep,"%02d",  ((int)ep[3]/12)*12); stat|=od_rtk_rtkcmn_repstr(rpath,"%hc",rep);
        sprintf(rep,"%04.0f",ep[0]);              stat|=od_rtk_rtkcmn_repstr(rpath,"%Y",rep);
        sprintf(rep,"%02.0f",fmod(ep[0],100.0));  stat|=od_rtk_rtkcmn_repstr(rpath,"%y",rep);
        sprintf(rep,"%02.0f",ep[1]);              stat|=od_rtk_rtkcmn_repstr(rpath,"%m",rep);
        sprintf(rep,"%02.0f",ep[2]);              stat|=od_rtk_rtkcmn_repstr(rpath,"%d",rep);
        sprintf(rep,"%02.0f",ep[3]);              stat|=od_rtk_rtkcmn_repstr(rpath,"%h",rep);
        sprintf(rep,"%02.0f",ep[4]);              stat|=od_rtk_rtkcmn_repstr(rpath,"%M",rep);
        sprintf(rep,"%02.0f",floor(ep[5]));       stat|=od_rtk_rtkcmn_repstr(rpath,"%S",rep);
        sprintf(rep,"%03d",  doy);                stat|=od_rtk_rtkcmn_repstr(rpath,"%n",rep);
        sprintf(rep,"%04d",  week);               stat|=od_rtk_rtkcmn_repstr(rpath,"%W",rep);
        sprintf(rep,"%d",    dow);                stat|=od_rtk_rtkcmn_repstr(rpath,"%D",rep);
        sprintf(rep,"%c",    'a'+(int)ep[3]);     stat|=od_rtk_rtkcmn_repstr(rpath,"%H",rep);
        sprintf(rep,"%02d",  ((int)ep[4]/15)*15); stat|=od_rtk_rtkcmn_repstr(rpath,"%t",rep);
    }
    else if (strstr(rpath,"%ha")||strstr(rpath,"%hb")||strstr(rpath,"%hc")||
             strstr(rpath,"%Y" )||strstr(rpath,"%y" )||strstr(rpath,"%m" )||
             strstr(rpath,"%d" )||strstr(rpath,"%h" )||strstr(rpath,"%M" )||
             strstr(rpath,"%S" )||strstr(rpath,"%n" )||strstr(rpath,"%W" )||
             strstr(rpath,"%D" )||strstr(rpath,"%H" )||strstr(rpath,"%t" )) {
        return -1; /* no valid time */
    }
    return stat;
}
/* replace keywords in file path and generate multiple paths -------------------
* replace keywords in file path with date, time, rover and base station id
* generate multiple keywords-replaced paths
* args   : char   *path     I   file path (see below)
*          char   *rpath[]  O   file paths in which keywords replaced
*          int    nmax      I   max number of output file paths
*          gtime_t ts       I   time start (gpst)
*          gtime_t te       I   time end   (gpst)
*          char   *rov      I   rover id string        ("": not replaced)
*          char   *base     I   base station id string ("": not replaced)
* return : number of replaced file paths
* notes  : see reppath() for replacements of keywords.
*          minimum interval of time replaced is 900s.
*-----------------------------------------------------------------------------*/
extern int reppaths(const char *path, char *rpath[], int nmax, gtime_t ts,
                    gtime_t te, const char *rov, const char *base)
{
    gtime_t time;
    double tow,tint=86400.0;
    int i,n=0,week;
    
    trace(3,"reppaths: path =%s nmax=%d rov=%s base=%s\n",path,nmax,rov,base);
    
    if (ts.time==0||te.time==0||timediff(ts,te)>0.0) return 0;
    
    if (strstr(path,"%S")||strstr(path,"%M")||strstr(path,"%t")) tint=900.0;
    else if (strstr(path,"%h")||strstr(path,"%H")) tint=3600.0;
    
    tow=time2gpst(ts,&week);
    time=gpst2time(week,floor(tow/tint)*tint);
    
    while (timediff(time,te)<=0.0&&n<nmax) {
        reppath(path,rpath[n],time,rov,base);
        if (n==0||strcmp(rpath[n],rpath[n-1])) n++;
        time=timeadd(time,tint);
    }
    for (i=0;i<n;i++) trace(3,"reppaths: rpath=%s\n",rpath[i]);
    return n;
}
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m)
*          double *rr       I   receiver position (ecef at reception) (m)
*          double *e        O   line-of-sight vector (ecef)
* return : geometric distance (m) (0>:error/no satellite position)
* notes  : distance includes sagnac effect correction
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
    double r;
    int i;
    
    if (norm(rs,3)<RE_WGS84) return -1.0;
    for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
    r=norm(e,3);
    for (i=0;i<3;i++) e[i]/=r;
    return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}
/* satellite azimuth/elevation angle -------------------------------------------
* compute satellite azimuth/elevation angle
* args   : double *pos      I   geodetic position {lat,lon,h} (rad,m)
*          double *e        I   receiver-to-satellilte unit vevtor (ecef)
*          double *azel     IO  azimuth/elevation {az,el} (rad) (NULL: no output)
*                               (0.0<=azel[0]<2*pi,-pi/2<=azel[1]<=pi/2)
* return : elevation angle (rad)
*-----------------------------------------------------------------------------*/
extern double satazel(const double *pos, const double *e, double *azel)
{
    double az=0.0,el=PI/2.0,enu[3];
    
    if (pos[2]>-RE_WGS84) {
        ecef2enu(pos,e,enu);
        az=dot(enu,enu,2)<1E-12?0.0:atan2(enu[0],enu[1]);
        if (az<0.0) az+=2*PI;
        el=asin(enu[2]);
    }
    if (azel) {azel[0]=az; azel[1]=el;}
    return el;
}
/* compute dops ----------------------------------------------------------------
* compute DOP (dilution of precision)
* args   : int    ns        I   number of satellites
*          double *azel     I   satellite azimuth/elevation angle (rad)
*          double elmin     I   elevation cutoff angle (rad)
*          double *dop      O   DOPs {GDOP,PDOP,HDOP,VDOP}
* return : none
* notes  : dop[0]-[3] return 0 in case of dop computation error
*-----------------------------------------------------------------------------*/
#define SQRT(x)     ((x)<0.0||(x)!=(x)?0.0:sqrt(x))

extern void dops(int ns, const double *azel, double elmin, double *dop)
{
    double H[4*MAXSAT],Q[16],cosel,sinel;
    int i,n;
    
    for (i=0;i<4;i++) dop[i]=0.0;
    for (i=n=0;i<ns&&i<MAXSAT;i++) {
        if (azel[1+i*2]<elmin||azel[1+i*2]<=0.0) continue;
        cosel=cos(azel[1+i*2]);
        sinel=sin(azel[1+i*2]);
        H[  4*n]=cosel*sin(azel[i*2]);
        H[1+4*n]=cosel*cos(azel[i*2]);
        H[2+4*n]=sinel;
        H[3+4*n++]=1.0;
    }
    if (n<4) return;
    
    matmul("NT",4,4,n,1.0,H,H,0.0,Q);
    if (!matinv(Q,4)) {
        dop[0]=SQRT(Q[0]+Q[5]+Q[10]+Q[15]); /* GDOP */
        dop[1]=SQRT(Q[0]+Q[5]+Q[10]);       /* PDOP */
        dop[2]=SQRT(Q[0]+Q[5]);             /* HDOP */
        dop[3]=SQRT(Q[10]);                 /* VDOP */
    }
}
/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
                       const double *azel)
{
    const double ion_default[]={ /* 2004/1/1 */
        0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
        0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
    };
    double tt,f,psi,phi,lam,amp,per,x;
    int week;
    
    if (pos[2]<-1E3||azel[1]<=0) return 0.0;
    if (norm(ion,8)<=0.0) ion=ion_default;
    
    /* earth centered angle (semi-circle) */
    psi=0.0137/(azel[1]/PI+0.11)-0.022;
    
    /* subionospheric latitude/longitude (semi-circle) */
    phi=pos[0]/PI+psi*cos(azel[0]);
    if      (phi> 0.416) phi= 0.416;
    else if (phi<-0.416) phi=-0.416;
    lam=pos[1]/PI+psi*sin(azel[0])/cos(phi*PI);
    
    /* geomagnetic latitude (semi-circle) */
    phi+=0.064*cos((lam-1.617)*PI);
    
    /* local time (s) */
    tt=43200.0*lam+time2gpst(t,&week);
    tt-=floor(tt/86400.0)*86400.0; /* 0<=tt<86400 */
    
    /* slant factor */
    f=1.0+16.0*pow(0.53-azel[1]/PI,3.0);
    
    /* ionospheric delay */
    amp=ion[0]+phi*(ion[1]+phi*(ion[2]+phi*ion[3]));
    per=ion[4]+phi*(ion[5]+phi*(ion[6]+phi*ion[7]));
    amp=amp<    0.0?    0.0:amp;
    per=per<72000.0?72000.0:per;
    x=2.0*PI*(tt-50400.0)/per;
    
    return CLIGHT*f*(fabs(x)<1.57?5E-9+amp*(1.0+x*x*(-0.5+x*x/24.0)):5E-9);
}
/* ionosphere mapping function -------------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
*-----------------------------------------------------------------------------*/
extern double ionmapf(const double *pos, const double *azel)
{
    if (pos[2]>=HION) return 1.0;
    return 1.0/cos(asin((RE_WGS84+pos[2])/(RE_WGS84+HION)*sin(PI/2.0-azel[1])));
}
/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double re        I   earth radius (km)
*          double hion      I   altitude of ionosphere (km)
*          double *posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
*-----------------------------------------------------------------------------*/
extern double ionppp(const double *pos, const double *azel, double re,
                     double hion, double *posp)
{
    double cosaz,rp,ap,sinap,tanap;
    
    rp=re/(re+hion)*cos(azel[1]);
    ap=PI/2.0-azel[1]-asin(rp);
    sinap=sin(ap);
    tanap=tan(ap);
    cosaz=cos(azel[0]);
    posp[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz);
    
    if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI/2.0-pos[0]))||
        (pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) {
        posp[1]=pos[1]+PI-asin(sinap*sin(azel[0])/cos(posp[0]));
    }
    else {
        posp[1]=pos[1]+asin(sinap*sin(azel[0])/cos(posp[0]));
    }
    return 1.0/sqrt(1.0-rp*rp);
}
/* troposphere model -----------------------------------------------------------
* compute tropospheric delay by standard atmosphere and saastamoinen model
* args   : gtime_t time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double humi      I   relative humidity
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double tropmodel(gtime_t time, const double *pos, const double *azel,
                        double humi)
{
    const double temp0=15.0; /* temparature at sea level */
    double hgt,pres,temp,e,z,trph,trpw;
    
    if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;
    
    /* standard atmosphere */
    hgt=pos[2]<0.0?0.0:pos[2];
    
    pres=1013.25*pow(1.0-2.2557E-5*hgt,5.2568);
    temp=temp0-6.5E-3*hgt+273.16;
    e=6.108*humi*exp((17.15*temp-4684.0)/(temp-38.45));
    
    /* saastamoninen model */
    z=PI/2.0-azel[1];
    trph=0.0022768*pres/(1.0-0.00266*cos(2.0*pos[0])-0.00028*hgt/1E3)/cos(z);
    trpw=0.002277*(1255.0/temp+0.05)*e/cos(z);
    return trph+trpw;
}
#ifndef IERS_MODEL

static double od_rtk_rtkcmn_interpc(const double coef[], double lat)
{
    int i=(int)(lat/15.0);
    if (i<1) return coef[0]; else if (i>4) return coef[4];
    return coef[i-1]*(1.0-lat/15.0+i)+coef[i]*(lat/15.0-i);
}
static double od_rtk_rtkcmn_mapf(double el, double a, double b, double c)
{
    double sinel=sin(el);
    return (1.0+a/(1.0+b/(1.0+c)))/(sinel+(a/(sinel+b/(sinel+c))));
}
static double od_rtk_rtkcmn_nmf(gtime_t time, const double pos[], const double azel[],
                  double *mapfw)
{
    /* ref [5] table 3 */
    /* hydro-ave-a,b,c, hydro-amp-a,b,c, wet-a,b,c at latitude 15,30,45,60,75 */
    const double coef[][5]={
        { 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3},
        { 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3},
        { 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3},
        
        { 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5},
        { 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5},
        { 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5},
        
        { 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4},
        { 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3},
        { 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2}
    };
    const double aht[]={ 2.53E-5, 5.49E-3, 1.14E-3}; /* height correction */
    
    double y,cosy,ah[3],aw[3],dm,el=azel[1],lat=pos[0]*R2D,hgt=pos[2];
    int i;
    
    if (el<=0.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }
    /* year from doy 28, added half a year for southern latitudes */
    y=(time2doy(time)-28.0)/365.25+(lat<0.0?0.5:0.0);
    
    cosy=cos(2.0*PI*y);
    lat=fabs(lat);
    
    for (i=0;i<3;i++) {
        ah[i]=od_rtk_rtkcmn_interpc(coef[i  ],lat)-od_rtk_rtkcmn_interpc(coef[i+3],lat)*cosy;
        aw[i]=od_rtk_rtkcmn_interpc(coef[i+6],lat);
    }
    /* ellipsoidal height is used instead of height above sea level */
    dm=(1.0/sin(el)-od_rtk_rtkcmn_mapf(el,aht[0],aht[1],aht[2]))*hgt/1E3;
    
    if (mapfw) *mapfw=od_rtk_rtkcmn_mapf(el,aw[0],aw[1],aw[2]);
    
    return od_rtk_rtkcmn_mapf(el,ah[0],ah[1],ah[2])+dm;
}
#endif /* !IERS_MODEL */

/* troposphere mapping function ------------------------------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
* note   : see ref [5] (NMF) and [9] (GMF)
*          original JGR paper of [5] has bugs in eq.(4) and (5). the corrected
*          paper is obtained from:
*          ftp://web.haystack.edu/pub/aen/nmf/NMF_JGR.pdf
*-----------------------------------------------------------------------------*/
extern double tropmapf(gtime_t time, const double pos[], const double azel[],
                       double *mapfw)
{
#ifdef IERS_MODEL
    const double ep[]={2000,1,1,12,0,0};
    double mjd,lat,lon,hgt,zd,gmfh,gmfw;
#endif
    trace(4,"tropmapf: pos=%10.6f %11.6f %6.1f azel=%5.1f %4.1f\n",
          pos[0]*R2D,pos[1]*R2D,pos[2],azel[0]*R2D,azel[1]*R2D);
    
    if (pos[2]<-1000.0||pos[2]>20000.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }
#ifdef IERS_MODEL
    mjd=51544.5+(timediff(time,epoch2time(ep)))/86400.0;
    lat=pos[0];
    lon=pos[1];
    hgt=pos[2]-geoidh(pos); /* height in m (mean sea level) */
    zd =PI/2.0-azel[1];
    
    /* call GMF */
    gmf_(&mjd,&lat,&lon,&hgt,&zd,&gmfh,&gmfw);
    
    if (mapfw) *mapfw=gmfw;
    return gmfh;
#else
    return od_rtk_rtkcmn_nmf(time,pos,azel,mapfw); /* NMF */
#endif
}
/* interpolate antenna phase center variation --------------------------------*/
static double od_rtk_rtkcmn_interpvar(double ang, const double *var)
{
    double a=ang/5.0; /* ang=0-90 */
    int i=(int)a;
    if (i<0) return var[0]; else if (i>=18) return var[18];
    return var[i]*(1.0-a+i)+var[i+1]*(a-i);
}
/* receiver antenna model ------------------------------------------------------
* compute antenna offset by antenna phase center parameters
* args   : pcv_t *pcv       I   antenna phase center parameters
*          double *del      I   antenna delta {e,n,u} (m)
*          double *azel     I   azimuth/elevation for receiver {az,el} (rad)
*          int     opt      I   option (0:only offset,1:offset+pcv)
*          double *dant     O   range offsets for each frequency (m)
* return : none
* notes  : current version does not support azimuth dependent terms
*-----------------------------------------------------------------------------*/
extern void antmodel(const pcv_t *pcv, const double *del, const double *azel,
                     int opt, double *dant)
{
    double e[3],off[3],cosel=cos(azel[1]);
    int i,j;
    
    trace(4,"antmodel: azel=%6.1f %4.1f opt=%d\n",azel[0]*R2D,azel[1]*R2D,opt);
    
    e[0]=sin(azel[0])*cosel;
    e[1]=cos(azel[0])*cosel;
    e[2]=sin(azel[1]);
    
    for (i=0;i<NFREQ;i++) {
        for (j=0;j<3;j++) off[j]=pcv->off[i][j]+del[j];
        
        dant[i]=-dot(off,e,3)+(opt?od_rtk_rtkcmn_interpvar(90.0-azel[1]*R2D,pcv->var[i]):0.0);
    }
    trace(5,"antmodel: dant=%6.3f %6.3f\n",dant[0],dant[1]);
}
/* satellite antenna model ------------------------------------------------------
* compute satellite antenna phase center parameters
* args   : pcv_t *pcv       I   antenna phase center parameters
*          double nadir     I   nadir angle for satellite (rad)
*          double *dant     O   range offsets for each frequency (m)
* return : none
*-----------------------------------------------------------------------------*/
extern void antmodel_s(const pcv_t *pcv, double nadir, double *dant)
{
    int i;
    
    trace(4,"antmodel_s: nadir=%6.1f\n",nadir*R2D);
    
    for (i=0;i<NFREQ;i++) {
        dant[i]=od_rtk_rtkcmn_interpvar(nadir*R2D*5.0,pcv->var[i]);
    }
    trace(5,"antmodel_s: dant=%6.3f %6.3f\n",dant[0],dant[1]);
}
/* sun and moon position in eci (ref [4] 5.1.1, 5.2.1) -----------------------*/
static void od_rtk_rtkcmn_sunmoonpos_eci(gtime_t tut, double *rsun, double *rmoon)
{
    const double ep2000[]={2000,1,1,12,0,0};
    double t,f[5],eps,Ms,ls,rs,lm,pm,rm,sine,cose,sinp,cosp,sinl,cosl;
    
    trace(4,"sunmoonpos_eci: tut=%s\n",time_str(tut,3));
    
    t=timediff(tut,epoch2time(ep2000))/86400.0/36525.0;
    
    /* astronomical arguments */
    od_rtk_rtkcmn_ast_args(t,f);
    
    /* obliquity of the ecliptic */
    eps=23.439291-0.0130042*t;
    sine=sin(eps*D2R); cose=cos(eps*D2R);
    
    /* sun position in eci */
    if (rsun) {
        Ms=357.5277233+35999.05034*t;
        ls=280.460+36000.770*t+1.914666471*sin(Ms*D2R)+0.019994643*sin(2.0*Ms*D2R);
        rs=AU*(1.000140612-0.016708617*cos(Ms*D2R)-0.000139589*cos(2.0*Ms*D2R));
        sinl=sin(ls*D2R); cosl=cos(ls*D2R);
        rsun[0]=rs*cosl;
        rsun[1]=rs*cose*sinl;
        rsun[2]=rs*sine*sinl;
        
        trace(5,"rsun =%.3f %.3f %.3f\n",rsun[0],rsun[1],rsun[2]);
    }
    /* moon position in eci */
    if (rmoon) {
        lm=218.32+481267.883*t+6.29*sin(f[0])-1.27*sin(f[0]-2.0*f[3])+
           0.66*sin(2.0*f[3])+0.21*sin(2.0*f[0])-0.19*sin(f[1])-0.11*sin(2.0*f[2]);
        pm=5.13*sin(f[2])+0.28*sin(f[0]+f[2])-0.28*sin(f[2]-f[0])-
           0.17*sin(f[2]-2.0*f[3]);
        rm=RE_WGS84/sin((0.9508+0.0518*cos(f[0])+0.0095*cos(f[0]-2.0*f[3])+
                   0.0078*cos(2.0*f[3])+0.0028*cos(2.0*f[0]))*D2R);
        sinl=sin(lm*D2R); cosl=cos(lm*D2R);
        sinp=sin(pm*D2R); cosp=cos(pm*D2R);
        rmoon[0]=rm*cosp*cosl;
        rmoon[1]=rm*(cose*cosp*sinl-sine*sinp);
        rmoon[2]=rm*(sine*cosp*sinl+cose*sinp);
        
        trace(5,"rmoon=%.3f %.3f %.3f\n",rmoon[0],rmoon[1],rmoon[2]);
    }
}
/* sun and moon position -------------------------------------------------------
* get sun and moon position in ecef
* args   : gtime_t tut      I   time in ut1
*          double *erpv     I   erp value {xp,yp,ut1_utc,lod} (rad,rad,s,s/d)
*          double *rsun     IO  sun position in ecef  (m) (NULL: not output)
*          double *rmoon    IO  moon position in ecef (m) (NULL: not output)
*          double *gmst     O   gmst (rad)
* return : none
*-----------------------------------------------------------------------------*/
extern void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
                       double *rmoon, double *gmst)
{
    gtime_t tut;
    double rs[3],rm[3],U[9],gmst_;
    
    trace(4,"sunmoonpos: tutc=%s\n",time_str(tutc,3));
    
    tut=timeadd(tutc,erpv[2]); /* utc -> ut1 */
    
    /* sun and moon position in eci */
    od_rtk_rtkcmn_sunmoonpos_eci(tut,rsun?rs:NULL,rmoon?rm:NULL);
    
    /* eci to ecef transformation matrix */
    eci2ecef(tutc,erpv,U,&gmst_);
    
    /* sun and moon postion in ecef */
    if (rsun ) matmul("NN",3,1,3,1.0,U,rs,0.0,rsun );
    if (rmoon) matmul("NN",3,1,3,1.0,U,rm,0.0,rmoon);
    if (gmst ) *gmst=gmst_;
}
/* uncompress file -------------------------------------------------------------
* uncompress (uncompress/unzip/uncompact hatanaka-compression/tar) file
* args   : char   *file     I   input file
*          char   *uncfile  O   uncompressed file
* return : status (-1:error,0:not compressed file,1:uncompress completed)
* note   : creates uncompressed file in tempolary directory
*          gzip, tar and crx2rnx commands have to be installed in commands path
*-----------------------------------------------------------------------------*/
extern int rtk_uncompress(const char *file, char *uncfile)
{
    int stat=0;
    char *p,cmd[64+2048]="",tmpfile[1024]="",buff[1024],*fname,*dir="";
    
    trace(3,"rtk_uncompress: file=%s\n",file);
    
    strcpy(tmpfile,file);
    if (!(p=strrchr(tmpfile,'.'))) return 0;
    
    /* uncompress by gzip */
    if (!strcmp(p,".z"  )||!strcmp(p,".Z"  )||
        !strcmp(p,".gz" )||!strcmp(p,".GZ" )||
        !strcmp(p,".zip")||!strcmp(p,".ZIP")) {
        
        strcpy(uncfile,tmpfile); uncfile[p-tmpfile]='\0';
        sprintf(cmd,"gzip -f -d -c \"%s\" > \"%s\"",tmpfile,uncfile);
        
        if (execcmd(cmd)) {
            remove(uncfile);
            return -1;
        }
        strcpy(tmpfile,uncfile);
        stat=1;
    }
    /* extract tar file */
    if ((p=strrchr(tmpfile,'.'))&&!strcmp(p,".tar")) {
        
        strcpy(uncfile,tmpfile); uncfile[p-tmpfile]='\0';
        strcpy(buff,tmpfile);
        fname=buff;
#ifdef WIN32
        if ((p=strrchr(buff,'\\'))) {
            *p='\0'; dir=fname; fname=p+1;
        }
        sprintf(cmd,"set PATH=%%CD%%;%%PATH%% & cd /D \"%s\" & tar -xf \"%s\"",
                dir,fname);
#else
        if ((p=strrchr(buff,'/'))) {
            *p='\0'; dir=fname; fname=p+1;
        }
        sprintf(cmd,"tar -C \"%s\" -xf \"%s\"",dir,tmpfile);
#endif
        if (execcmd(cmd)) {
            if (stat) remove(tmpfile);
            return -1;
        }
        if (stat) remove(tmpfile);
        stat=1;
    }
    /* extract hatanaka-compressed file by cnx2rnx */
    else if ((p=strrchr(tmpfile,'.'))&&
             ((strlen(p)>3&&(*(p+3)=='d'||*(p+3)=='D'))||
              !strcmp(p,".crx")||!strcmp(p,".CRX"))) {
        
        strcpy(uncfile,tmpfile);
        uncfile[p-tmpfile+3]=*(p+3)=='D'?'O':'o';
        sprintf(cmd,"crx2rnx < \"%s\" > \"%s\"",tmpfile,uncfile);
        
        if (execcmd(cmd)) {
            remove(uncfile);
            if (stat) remove(tmpfile);
            return -1;
        }
        if (stat) remove(tmpfile);
        stat=1;
    }
    trace(3,"rtk_uncompress: stat=%d\n",stat);
    return stat;
}
/* dummy application functions for shared library ----------------------------*/
#ifdef WIN_DLL
extern int showmsg(const char *format,...) {return 0;}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}
#endif



#pragma pop_macro("dgetrs_")
#pragma pop_macro("dgetri_")
#pragma pop_macro("dgetrf_")
#pragma pop_macro("dgemm_")
#pragma pop_macro("_POSIX_C_SOURCE")
#pragma pop_macro("SQRT")
#pragma pop_macro("SQR")
#pragma pop_macro("Rz")
#pragma pop_macro("Ry")
#pragma pop_macro("Rx")
#pragma pop_macro("POLYCRC32")
#pragma pop_macro("POLYCRC24Q")
#pragma pop_macro("MAX_VAR_EPH")
#pragma pop_macro("LAPACK")


/* ===== Embedded rtcm.c ===== */
#pragma push_macro("RTCM2PREAMB")
#undef RTCM2PREAMB
#pragma push_macro("RTCM3PREAMB")
#undef RTCM3PREAMB

/*------------------------------------------------------------------------------
* rtcm.c : rtcm functions
*
*          Copyright (C) 2009-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1]  RTCM Recommended Standards for Differential GNSS (Global Navigation
*          Satellite Systems) Service version 2.3, August 20, 2001
*     [7]  RTCM Standard 10403.1 - Amendment 5, Differential GNSS (Global
*          Navigation Satellite Systems) Services - version 3, July 1, 2011
*     [10] RTCM Paper 059-2011-SC104-635 (draft Galileo and QZSS ssr messages)
*     [15] RTCM Standard 10403.2, Differential GNSS (Global Navigation Satellite
*          Systems) Services - version 3, with amendment 1/2, November 7, 2013
*     [16] Proposal of new RTCM SSR Messages (ssr_1_gal_qzss_sbas_dbs_v05)
*          2014/04/17
*     [17] RTCM Standard 10403.3, Differential GNSS (Global Navigation Satellite
*          Systems) Services - version 3, with amendment 1, April 28, 2020
*     [18] IGS State Space Representation (SSR) Format version 1.00, October 5,
*          2020
*
* version : $Revision:$ $Date:$
* history : 2009/04/10 1.0  new
*           2009/06/29 1.1  support type 1009-1012 to get synchronous-gnss-flag
*           2009/12/04 1.2  support type 1010,1012,1020
*           2010/07/15 1.3  support type 1057-1068 for ssr corrections
*                           support type 1007,1008,1033 for antenna info
*           2010/09/08 1.4  fix problem of ephemeris and ssr sequence upset
*                           (2.4.0_p8)
*           2012/05/11 1.5  comply with RTCM 3 final SSR format (RTCM 3
*                           Amendment 5) (ref [7]) (2.4.1_p6)
*           2012/05/14 1.6  separate rtcm2.c, rtcm3.c
*                           add options to select used codes for msm
*           2013/04/27 1.7  comply with rtcm 3.2 with amendment 1/2 (ref[15])
*           2013/12/06 1.8  support SBAS/BeiDou SSR messages (ref[16])
*           2018/01/29 1.9  support RTCM 3.3 (ref[17])
*                           crc24q() -> rtk_crc24q()
*           2018/10/10 1.10 fix bug on initializing rtcm struct
*                           add rtcm option -GALINAV, -GALFNAV
*           2018/11/05 1.11 add notes for api gen_rtcm3()
*           2020/11/30 1.12 modify API gen_rtcm3()
*                           support NavIC/IRNSS MSM and ephemeris (ref [17])
*                           allocate double size of ephemeris buffer to support
*                            multiple ephemeris sets in init_rtcm()
*                           delete references [2]-[6],[8],[9],[11]-[14]
*                           update reference [17]
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/

/* function prototypes -------------------------------------------------------*/
extern int decode_rtcm2(rtcm_t *rtcm);
extern int decode_rtcm3(rtcm_t *rtcm);
extern int encode_rtcm3(rtcm_t *rtcm, int type, int subtype, int sync);

/* constants -----------------------------------------------------------------*/

#define RTCM2PREAMB 0x66        /* rtcm ver.2 frame preamble */
#define RTCM3PREAMB 0xD3        /* rtcm ver.3 frame preamble */

/* initialize rtcm control -----------------------------------------------------
* initialize rtcm control struct and reallocate memory for observation and
* ephemeris buffer in rtcm control struct
* args   : rtcm_t *raw      IO  rtcm control struct
* return : status (1:ok,0:memory allocation error)
*-----------------------------------------------------------------------------*/
extern int init_rtcm(rtcm_t *rtcm)
{
    gtime_t time0={0};
    obsd_t data0={{0}};
    eph_t  eph0 ={0,-1,-1};
    geph_t geph0={0,-1};
    ssr_t ssr0={{{0}}};
    int i,j;
    
    trace(3,"init_rtcm:\n");
    
    rtcm->staid=rtcm->stah=rtcm->seqno=rtcm->outtype=0;
    rtcm->time=rtcm->time_s=time0;
    rtcm->sta.name[0]=rtcm->sta.marker[0]='\0';
    rtcm->sta.antdes[0]=rtcm->sta.antsno[0]='\0';
    rtcm->sta.rectype[0]=rtcm->sta.recver[0]=rtcm->sta.recsno[0]='\0';
    rtcm->sta.antsetup=rtcm->sta.itrf=rtcm->sta.deltype=0;
    for (i=0;i<3;i++) {
        rtcm->sta.pos[i]=rtcm->sta.del[i]=0.0;
    }
    rtcm->sta.hgt=0.0;
    rtcm->dgps=NULL;
    for (i=0;i<MAXSAT;i++) {
        rtcm->ssr[i]=ssr0;
    }
    rtcm->msg[0]=rtcm->msgtype[0]=rtcm->opt[0]='\0';
    for (i=0;i<6;i++) rtcm->msmtype[i][0]='\0';
    rtcm->obsflag=rtcm->ephsat=0;
    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ+NEXOBS;j++) {
        rtcm->cp[i][j]=0.0;
        rtcm->lock[i][j]=rtcm->loss[i][j]=0;
        rtcm->lltime[i][j]=time0;
    }
    rtcm->nbyte=rtcm->nbit=rtcm->len=0;
    rtcm->word=0;
    for (i=0;i<100;i++) rtcm->nmsg2[i]=0;
    for (i=0;i<400;i++) rtcm->nmsg3[i]=0;
    
    rtcm->obs.data=NULL;
    rtcm->nav.eph =NULL;
    rtcm->nav.geph=NULL;
    
    /* reallocate memory for observation and ephemeris buffer */
    if (!(rtcm->obs.data=(obsd_t *)malloc(sizeof(obsd_t)*MAXOBS))||
        !(rtcm->nav.eph =(eph_t  *)malloc(sizeof(eph_t )*MAXSAT*2))||
        !(rtcm->nav.geph=(geph_t *)malloc(sizeof(geph_t)*MAXPRNGLO))) {
        free_rtcm(rtcm);
        return 0;
    }
    rtcm->obs.n=0;
    rtcm->nav.n=MAXSAT*2;
    rtcm->nav.ng=MAXPRNGLO;
    for (i=0;i<MAXOBS   ;i++) rtcm->obs.data[i]=data0;
    for (i=0;i<MAXSAT*2 ;i++) rtcm->nav.eph [i]=eph0;
    for (i=0;i<MAXPRNGLO;i++) rtcm->nav.geph[i]=geph0;
    return 1;
}
/* free rtcm control ----------------------------------------------------------
* free observation and ephemeris buffer in rtcm control struct
* args   : rtcm_t *raw      IO  rtcm control struct
* return : none
*-----------------------------------------------------------------------------*/
extern void free_rtcm(rtcm_t *rtcm)
{
    trace(3,"free_rtcm:\n");
    
    /* free memory for observation and ephemeris buffer */
    free(rtcm->obs.data); rtcm->obs.data=NULL; rtcm->obs.n=0;
    free(rtcm->nav.eph ); rtcm->nav.eph =NULL; rtcm->nav.n=0;
    free(rtcm->nav.geph); rtcm->nav.geph=NULL; rtcm->nav.ng=0;
}
/* input RTCM 2 message from stream --------------------------------------------
* fetch next RTCM 2 message and input a message from byte stream
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          uint8_t data     I   stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 5: input station pos/ant parameters,
*                  6: input time parameter, 7: input dgps corrections,
*                  9: input special message)
* notes  : before firstly calling the function, time in rtcm control struct has
*          to be set to the approximate time within 1/2 hour in order to resolve
*          ambiguity of time in rtcm messages.
*          supported msgs RTCM ver.2: 1,3,9,14,16,17,18,19,22
*          refer [1] for RTCM ver.2
*-----------------------------------------------------------------------------*/
extern int input_rtcm2(rtcm_t *rtcm, uint8_t data)
{
    uint8_t preamb;
    int i;
    
    trace(5,"input_rtcm2: data=%02x\n",data);
    
    if ((data&0xC0)!=0x40) return 0; /* ignore if upper 2bit != 01 */
    
    for (i=0;i<6;i++,data>>=1) { /* decode 6-of-8 form */
        rtcm->word=(rtcm->word<<1)+(data&1);
        
        /* synchronize frame */
        if (rtcm->nbyte==0) {
            preamb=(uint8_t)(rtcm->word>>22);
            if (rtcm->word&0x40000000) preamb^=0xFF; /* decode preamble */
            if (preamb!=RTCM2PREAMB) continue;
            
            /* check parity */
            if (!decode_word(rtcm->word,rtcm->buff)) continue;
            rtcm->nbyte=3; rtcm->nbit=0;
            continue;
        }
        if (++rtcm->nbit<30) continue; else rtcm->nbit=0;
        
        /* check parity */
        if (!decode_word(rtcm->word,rtcm->buff+rtcm->nbyte)) {
            trace(2,"rtcm2 partity error: i=%d word=%08x\n",i,rtcm->word);
            rtcm->nbyte=0; rtcm->word&=0x3;
            continue;
        }
        rtcm->nbyte+=3;
        if (rtcm->nbyte==6) rtcm->len=(rtcm->buff[5]>>3)*3+6;
        if (rtcm->nbyte<rtcm->len) continue;
        rtcm->nbyte=0; rtcm->word&=0x3;
        
        /* decode rtcm2 message */
        return decode_rtcm2(rtcm);
    }
    return 0;
}
/* input RTCM 3 message from stream --------------------------------------------
* fetch next RTCM 3 message and input a message from byte stream
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          uint8_t data     I   stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 5: input station pos/ant parameters,
*                  10: input ssr messages)
* notes  : before firstly calling the function, time in rtcm control struct has
*          to be set to the approximate time within 1/2 week in order to resolve
*          ambiguity of time in rtcm messages.
*          
*          to specify input options, set rtcm->opt to the following option
*          strings separated by spaces.
*
*          -EPHALL  : input all ephemerides (default: only new)
*          -STA=nnn : input only message with STAID=nnn (default: all)
*          -GLss    : select signal ss for GPS MSM (ss=1C,1P,...)
*          -RLss    : select signal ss for GLO MSM (ss=1C,1P,...)
*          -ELss    : select signal ss for GAL MSM (ss=1C,1B,...)
*          -JLss    : select signal ss for QZS MSM (ss=1C,2C,...)
*          -CLss    : select signal ss for BDS MSM (ss=2I,7I,...)
*          -ILss    : select signal ss for IRN MSM (ss=5A,9A,...)
*          -GALINAV : select I/NAV for Galileo ephemeris (default: all)
*          -GALFNAV : select F/NAV for Galileo ephemeris (default: all)
*
*          supported RTCM 3 messages (ref [7][10][15][16][17][18])
*
*            TYPE       :  GPS   GLONASS Galileo  QZSS     BDS    SBAS    NavIC
*         ----------------------------------------------------------------------
*          OBS COMP L1  : 1001~   1009~     -       -       -       -       -
*              FULL L1  : 1002    1010      -       -       -       -       -
*              COMP L1L2: 1003~   1011~     -       -       -       -       -
*              FULL L1L2: 1004    1012      -       -       -       -       -
*
*          NAV          : 1019    1020    1045**  1044    1042      -     1041
*                           -       -     1046**    -       63*     -       -
*
*          MSM 1        : 1071~   1081~   1091~   1111~   1121~   1101~   1131~
*              2        : 1072~   1082~   1092~   1112~   1122~   1102~   1132~
*              3        : 1073~   1083~   1093~   1113~   1123~   1103~   1133~
*              4        : 1074    1084    1094    1114    1124    1104    1134
*              5        : 1075    1085    1095    1115    1125    1105    1135 
*              6        : 1076    1086    1096    1116    1126    1106    1136 
*              7        : 1077    1087    1097    1117    1127    1107    1137 
*
*          SSR ORBIT    : 1057    1063    1240*   1246*   1258*     -       -
*              CLOCK    : 1058    1064    1241*   1247*   1259*     -       -
*              CODE BIAS: 1059    1065    1242*   1248*   1260*     -       -
*              OBT/CLK  : 1060    1066    1243*   1249*   1261*     -       -
*              URA      : 1061    1067    1244*   1250*   1262*     -       -
*              HR-CLOCK : 1062    1068    1245*   1251*   1263*     -       -
*              PHAS BIAS:   11*     -       12*     13*     14*     -       -
*
*          ANT/RCV INFO : 1007    1008    1033
*          STA POSITION : 1005    1006
*
*          PROPRIETARY  : 4076 (IGS)
*         ----------------------------------------------------------------------
*                            (* draft, ** 1045:F/NAV,1046:I/NAV, ~ only encode)
*
*          for MSM observation data with multiple signals for a frequency,
*          a signal is selected according to internal priority. to select
*          a specified signal, use the input options.
*
*          RTCM 3 message format:
*            +----------+--------+-----------+--------------------+----------+
*            | preamble | 000000 |  length   |    data message    |  parity  |
*            +----------+--------+-----------+--------------------+----------+
*            |<-- 8 --->|<- 6 -->|<-- 10 --->|<--- length x 8 --->|<-- 24 -->|
*            
*-----------------------------------------------------------------------------*/
extern int input_rtcm3(rtcm_t *rtcm, uint8_t data)
{
    trace(5,"input_rtcm3: data=%02x\n",data);
    
    /* synchronize frame */
    if (rtcm->nbyte==0) {
        if (data!=RTCM3PREAMB) return 0;
        rtcm->buff[rtcm->nbyte++]=data;
        return 0;
    }
    rtcm->buff[rtcm->nbyte++]=data;
    
    if (rtcm->nbyte==3) {
        rtcm->len=getbitu(rtcm->buff,14,10)+3; /* length without parity */
    }
    if (rtcm->nbyte<3||rtcm->nbyte<rtcm->len+3) return 0;
    rtcm->nbyte=0;
    
    /* check parity */
    if (rtk_crc24q(rtcm->buff,rtcm->len)!=getbitu(rtcm->buff,rtcm->len*8,24)) {
        trace(2,"rtcm3 parity error: len=%d\n",rtcm->len);
        return 0;
    }
    /* decode rtcm3 message */
    return decode_rtcm3(rtcm);
}
/* input RTCM 2 message from file ----------------------------------------------
* fetch next RTCM 2 message and input a messsage from file
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          FILE  *fp        I   file pointer
* return : status (-2: end of file, -1...10: same as above)
* notes  : same as above
*-----------------------------------------------------------------------------*/
extern int input_rtcm2f(rtcm_t *rtcm, FILE *fp)
{
    int i,data=0,ret;
    
    trace(4,"input_rtcm2f: data=%02x\n",data);
    
    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_rtcm2(rtcm,(uint8_t)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}
/* input RTCM 3 message from file ----------------------------------------------
* fetch next RTCM 3 message and input a messsage from file
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          FILE  *fp        I   file pointer
* return : status (-2: end of file, -1...10: same as above)
* notes  : same as above
*-----------------------------------------------------------------------------*/
extern int input_rtcm3f(rtcm_t *rtcm, FILE *fp)
{
    int i,data=0,ret;
    
    trace(4,"input_rtcm3f: data=%02x\n",data);
    
    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_rtcm3(rtcm,(uint8_t)data))) return ret;
    }
    return 0; /* return at every 4k bytes */
}
/* generate RTCM 2 message -----------------------------------------------------
* generate RTCM 2 message
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          int    type      I   message type
*          int    sync      I   sync flag (1:another message follows)
* return : status (1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int gen_rtcm2(rtcm_t *rtcm, int type, int sync)
{
    trace(4,"gen_rtcm2: type=%d sync=%d\n",type,sync);
    
    rtcm->nbit=rtcm->len=rtcm->nbyte=0;
    
    /* not yet implemented */
    
    return 0;
}
/* generate RTCM 3 message -----------------------------------------------------
* generate RTCM 3 message
* args   : rtcm_t *rtcm     IO  rtcm control struct
*          int    type      I   message type
*          int    subtype   I   message subtype
*          int    sync      I   sync flag (1:another message follows)
* return : status (1:ok,0:error)
* notes  : For rtcm 3 msm, the {nsat} x {nsig} in rtcm->obs should not exceed
*          64. If {nsat} x {nsig} of the input obs data exceeds 64, separate
*          them to multiple ones and call gen_rtcm3() multiple times as user
*          responsibility.
*          ({nsat} = number of valid satellites, {nsig} = number of signals in
*          the obs data) 
*-----------------------------------------------------------------------------*/
extern int gen_rtcm3(rtcm_t *rtcm, int type, int subtype, int sync)
{
    uint32_t crc;
    int i=0;
    
    trace(4,"gen_rtcm3: type=%d subtype=%d sync=%d\n",type,subtype,sync);
    
    rtcm->nbit=rtcm->len=rtcm->nbyte=0;
    
    /* set preamble and reserved */
    setbitu(rtcm->buff,i, 8,RTCM3PREAMB); i+= 8;
    setbitu(rtcm->buff,i, 6,0          ); i+= 6;
    setbitu(rtcm->buff,i,10,0          ); i+=10;
    
    /* encode rtcm 3 message body */
    if (!encode_rtcm3(rtcm,type,subtype,sync)) return 0;
    
    /* padding to align 8 bit boundary */
    for (i=rtcm->nbit;i%8;i++) {
        setbitu(rtcm->buff,i,1,0);
    }
    /* message length (header+data) (bytes) */
    if ((rtcm->len=i/8)>=3+1024) {
        trace(2,"generate rtcm 3 message length error len=%d\n",rtcm->len-3);
        rtcm->nbit=rtcm->len=0;
        return 0;
    }
    /* message length without header and parity */
    setbitu(rtcm->buff,14,10,rtcm->len-3);
    
    /* crc-24q */
    crc=rtk_crc24q(rtcm->buff,rtcm->len);
    setbitu(rtcm->buff,i,24,crc);
    
    /* length total (bytes) */
    rtcm->nbyte=rtcm->len+3;
    
    return 1;
}


#pragma pop_macro("RTCM3PREAMB")
#pragma pop_macro("RTCM2PREAMB")


/* ===== Embedded rtcm2.c ===== */

/*------------------------------------------------------------------------------
* rtcm2.c : rtcm ver.2 message functions
*
*          Copyright (C) 2009-2014 by T.TAKASU, All rights reserved.
*
* references :
*     see rtcm.c
*
* version : $Revision:$ $Date:$
* history : 2011/11/28 1.0  separated from rtcm.c
*           2014/10/21 1.1  fix problem on week rollover in rtcm 2 type 14
*-----------------------------------------------------------------------------*/

/* adjust hourly rollover of rtcm 2 time -------------------------------------*/
static void od_rtk_rtcm2_adjhour(rtcm_t *rtcm, double zcnt)
{
    double tow,hour,sec;
    int week;
    
    /* if no time, get cpu time */
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tow=time2gpst(rtcm->time,&week);
    hour=floor(tow/3600.0);
    sec=tow-hour*3600.0;
    if      (zcnt<sec-1800.0) zcnt+=3600.0;
    else if (zcnt>sec+1800.0) zcnt-=3600.0;
    rtcm->time=gpst2time(week,hour*3600+zcnt);
}
/* get observation data index ------------------------------------------------*/
static int od_rtk_rtcm2_obsindex(obs_t *obs, gtime_t time, int sat)
{
    int i,j;
    
    for (i=0;i<obs->n;i++) {
        if (obs->data[i].sat==sat) return i; /* field already exists */
    }
    if (i>=MAXOBS) return -1; /* overflow */
    
    /* add new field */
    obs->data[i].time=time;
    obs->data[i].sat=sat;
    for (j=0;j<NFREQ;j++) {
        obs->data[i].L[j]=obs->data[i].P[j]=0.0;
        obs->data[i].D[j]=0.0;
        obs->data[i].SNR[j]=obs->data[i].LLI[j]=obs->data[i].code[j]=0;
    }
    obs->n++;
    return i;
}
/* decode type 1/9: differential gps correction/partial correction set -------*/
static int od_rtk_rtcm2_decode_type1(rtcm_t *rtcm)
{
    int i=48,fact,udre,prn,sat,iod;
    double prc,rrc;
    
    trace(4,"decode_type1: len=%d\n",rtcm->len);
    
    while (i+40<=rtcm->len*8) {
        fact=getbitu(rtcm->buff,i, 1); i+= 1;
        udre=getbitu(rtcm->buff,i, 2); i+= 2;
        prn =getbitu(rtcm->buff,i, 5); i+= 5;
        prc =getbits(rtcm->buff,i,16); i+=16;
        rrc =getbits(rtcm->buff,i, 8); i+= 8;
        iod =getbits(rtcm->buff,i, 8); i+= 8;
        if (prn==0) prn=32;
        if (prc==0x80000000||rrc==0xFFFF8000) {
            trace(2,"rtcm2 1 prc/rrc indicates satellite problem: prn=%d\n",prn);
            continue;
        }
        if (rtcm->dgps) {
            sat=satno(SYS_GPS,prn);
            rtcm->dgps[sat-1].t0=rtcm->time;
            rtcm->dgps[sat-1].prc=prc*(fact?0.32:0.02);
            rtcm->dgps[sat-1].rrc=rrc*(fact?0.032:0.002);
            rtcm->dgps[sat-1].iod=iod;
            rtcm->dgps[sat-1].udre=udre;
        }
    }
    return 7;
}
/* decode type 3: reference station parameter --------------------------------*/
static int od_rtk_rtcm2_decode_type3(rtcm_t *rtcm)
{
    int i=48;
    
    trace(4,"decode_type3: len=%d\n",rtcm->len);
    
    if (i+96<=rtcm->len*8) {
        rtcm->sta.pos[0]=getbits(rtcm->buff,i,32)*0.01; i+=32;
        rtcm->sta.pos[1]=getbits(rtcm->buff,i,32)*0.01; i+=32;
        rtcm->sta.pos[2]=getbits(rtcm->buff,i,32)*0.01;
    }
    else {
        trace(2,"rtcm2 3 length error: len=%d\n",rtcm->len);
        return -1;
    }
    return 5;
}
/* decode type 14: gps time of week ------------------------------------------*/
static int od_rtk_rtcm2_decode_type14(rtcm_t *rtcm)
{
    double zcnt;
    int i=48,week,hour,leaps;
    
    trace(4,"decode_type14: len=%d\n",rtcm->len);
    
    zcnt=getbitu(rtcm->buff,24,13);
    if (i+24<=rtcm->len*8) {
        week =getbitu(rtcm->buff,i,10); i+=10;
        hour =getbitu(rtcm->buff,i, 8); i+= 8;
        leaps=getbitu(rtcm->buff,i, 6);
    }
    else {
        trace(2,"rtcm2 14 length error: len=%d\n",rtcm->len);
        return -1;
    }
    week=adjgpsweek(week);
    rtcm->time=gpst2time(week,hour*3600.0+zcnt*0.6);
    rtcm->nav.utc_gps[4]=leaps;
    return 6;
}
/* decode type 16: gps special message ---------------------------------------*/
static int od_rtk_rtcm2_decode_type16(rtcm_t *rtcm)
{
    int i=48,n=0;
    
    trace(4,"decode_type16: len=%d\n",rtcm->len);
    
    while (i+8<=rtcm->len*8&&n<90) {
        rtcm->msg[n++]=getbitu(rtcm->buff,i,8); i+=8;
    }
    rtcm->msg[n]='\0';
    
    trace(3,"rtcm2 16 message: %s\n",rtcm->msg);
    return 9;
}
/* decode type 17: gps ephemerides -------------------------------------------*/
static int od_rtk_rtcm2_decode_type17(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA;
    int i=48,week,prn,sat;
    
    trace(4,"decode_type17: len=%d\n",rtcm->len);
    
    if (i+480<=rtcm->len*8) {
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16);              i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16);              i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.iodc  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        prn       =getbitu(rtcm->buff,i, 5);              i+= 5+3;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.code  =getbitu(rtcm->buff,i, 2);              i+= 2;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.svh   =getbitu(rtcm->buff,i, 6);              i+= 6;
        eph.flag  =getbitu(rtcm->buff,i, 1);
    }
    else {
        trace(2,"rtcm2 17 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (prn==0) prn=32;
    sat=satno(SYS_GPS,prn);
    eph.sat=sat;
    eph.week=adjgpsweek(week);
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    rtcm->nav.eph[sat-1]=eph;
    rtcm->ephset=0;
    rtcm->ephsat=sat;
    return 2;
}
/* decode type 18: rtk uncorrected carrier-phase -----------------------------*/
static int od_rtk_rtcm2_decode_type18(rtcm_t *rtcm)
{
    gtime_t time;
    double usec,cp,tt;
    int i=48,index,freq,sync=1,code,sys,prn,sat,loss;
    
    trace(4,"decode_type18: len=%d\n",rtcm->len);
    
    if (i+24<=rtcm->len*8) {
        freq=getbitu(rtcm->buff,i, 2); i+= 2+2;
        usec=getbitu(rtcm->buff,i,20); i+=20;
    }
    else {
        trace(2,"rtcm2 18 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (freq&0x1) {
        trace(2,"rtcm2 18 not supported frequency: freq=%d\n",freq);
        return -1;
    }
    freq>>=1;
    
    while (i+48<=rtcm->len*8&&rtcm->obs.n<MAXOBS) {
        sync=getbitu(rtcm->buff,i, 1); i+= 1;
        code=getbitu(rtcm->buff,i, 1); i+= 1;
        sys =getbitu(rtcm->buff,i, 1); i+= 1;
        prn =getbitu(rtcm->buff,i, 5); i+= 5+3;
        loss=getbitu(rtcm->buff,i, 5); i+= 5;
        cp  =getbits(rtcm->buff,i,32); i+=32;
        if (prn==0) prn=32;
        if (!(sat=satno(sys?SYS_GLO:SYS_GPS,prn))) {
            trace(2,"rtcm2 18 satellite number error: sys=%d prn=%d\n",sys,prn);
            continue;
        }
        time=timeadd(rtcm->time,usec*1E-6);
        if (sys) time=utc2gpst(time); /* convert glonass time to gpst */
        
        tt=timediff(rtcm->obs.data[0].time,time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=od_rtk_rtcm2_obsindex(&rtcm->obs,time,sat))>=0) {
            rtcm->obs.data[index].L[freq]=-cp/256.0;
            rtcm->obs.data[index].LLI[freq]=rtcm->loss[sat-1][freq]!=loss;
            rtcm->obs.data[index].code[freq]=
                !freq?(code?CODE_L1P:CODE_L1C):(code?CODE_L2P:CODE_L2C);
            rtcm->loss[sat-1][freq]=loss;
        }
    }
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 19: rtk uncorrected pseudorange -------------------------------*/
static int od_rtk_rtcm2_decode_type19(rtcm_t *rtcm)
{
    gtime_t time;
    double usec,pr,tt;
    int i=48,index,freq,sync=1,code,sys,prn,sat;
    
    trace(4,"decode_type19: len=%d\n",rtcm->len);
    
    if (i+24<=rtcm->len*8) {
        freq=getbitu(rtcm->buff,i, 2); i+= 2+2;
        usec=getbitu(rtcm->buff,i,20); i+=20;
    }
    else {
        trace(2,"rtcm2 19 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (freq&0x1) {
        trace(2,"rtcm2 19 not supported frequency: freq=%d\n",freq);
        return -1;
    }
    freq>>=1;
    
    while (i+48<=rtcm->len*8&&rtcm->obs.n<MAXOBS) {
        sync=getbitu(rtcm->buff,i, 1); i+= 1;
        code=getbitu(rtcm->buff,i, 1); i+= 1;
        sys =getbitu(rtcm->buff,i, 1); i+= 1;
        prn =getbitu(rtcm->buff,i, 5); i+= 5+8;
        pr  =getbitu(rtcm->buff,i,32); i+=32;
        if (prn==0) prn=32;
        if (!(sat=satno(sys?SYS_GLO:SYS_GPS,prn))) {
            trace(2,"rtcm2 19 satellite number error: sys=%d prn=%d\n",sys,prn);
            continue;
        }
        time=timeadd(rtcm->time,usec*1E-6);
        if (sys) time=utc2gpst(time); /* convert glonass time to gpst */
        
        tt=timediff(rtcm->obs.data[0].time,time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=od_rtk_rtcm2_obsindex(&rtcm->obs,time,sat))>=0) {
            rtcm->obs.data[index].P[freq]=pr*0.02;
            rtcm->obs.data[index].code[freq]=
                !freq?(code?CODE_L1P:CODE_L1C):(code?CODE_L2P:CODE_L2C);
        }
    }
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 22: extended reference station parameter ----------------------*/
static int od_rtk_rtcm2_decode_type22(rtcm_t *rtcm)
{
    double del[2][3]={{0}},hgt=0.0;
    int i=48,j,noh;
    
    trace(4,"decode_type22: len=%d\n",rtcm->len);
    
    if (i+24<=rtcm->len*8) {
        del[0][0]=getbits(rtcm->buff,i,8)/25600.0; i+=8;
        del[0][1]=getbits(rtcm->buff,i,8)/25600.0; i+=8;
        del[0][2]=getbits(rtcm->buff,i,8)/25600.0; i+=8;
    }
    else {
        trace(2,"rtcm2 22 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (i+24<=rtcm->len*8) {
        i+=5; noh=getbits(rtcm->buff,i,1); i+=1;
        hgt=noh?0.0:getbitu(rtcm->buff,i,18)/25600.0;
        i+=18;
    }
    if (i+24<=rtcm->len*8) {
        del[1][0]=getbits(rtcm->buff,i,8)/1600.0; i+=8;
        del[1][1]=getbits(rtcm->buff,i,8)/1600.0; i+=8;
        del[1][2]=getbits(rtcm->buff,i,8)/1600.0;
    }
    rtcm->sta.deltype=1; /* xyz */
    for (j=0;j<3;j++) rtcm->sta.del[j]=del[0][j];
    rtcm->sta.hgt=hgt;
    return 5;
}
/* decode type 23: antenna type definition record ----------------------------*/
static int od_rtk_rtcm2_decode_type23(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 24: antenna reference point (arp) -----------------------------*/
static int od_rtk_rtcm2_decode_type24(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 31: differential glonass correction ---------------------------*/
static int od_rtk_rtcm2_decode_type31(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 32: differential glonass reference station parameters ---------*/
static int od_rtk_rtcm2_decode_type32(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 34: glonass partial differential correction set ---------------*/
static int od_rtk_rtcm2_decode_type34(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 36: glonass special message -----------------------------------*/
static int od_rtk_rtcm2_decode_type36(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 37: gnss system time offset -----------------------------------*/
static int od_rtk_rtcm2_decode_type37(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 59: proprietary message ---------------------------------------*/
static int od_rtk_rtcm2_decode_type59(rtcm_t *rtcm)
{
    return 0;
}
/* decode rtcm ver.2 message -------------------------------------------------*/
extern int decode_rtcm2(rtcm_t *rtcm)
{
    double zcnt;
    int staid,seqno,stah,ret=0,type=getbitu(rtcm->buff,8,6);
    
    trace(3,"decode_rtcm2: type=%2d len=%3d\n",type,rtcm->len);
    
    if ((zcnt=getbitu(rtcm->buff,24,13)*0.6)>=3600.0) {
        trace(2,"rtcm2 modified z-count error: zcnt=%.1f\n",zcnt);
        return -1;
    }
    od_rtk_rtcm2_adjhour(rtcm,zcnt);
    staid=getbitu(rtcm->buff,14,10);
    seqno=getbitu(rtcm->buff,37, 3);
    stah =getbitu(rtcm->buff,45, 3);
    if (seqno-rtcm->seqno!=1&&seqno-rtcm->seqno!=-7) {
        trace(2,"rtcm2 message outage: seqno=%d->%d\n",rtcm->seqno,seqno);
    }
    rtcm->seqno=seqno;
    rtcm->stah =stah;
    
    if (rtcm->outtype) {
        sprintf(rtcm->msgtype,"RTCM %2d (%4d) zcnt=%7.1f staid=%3d seqno=%d",
                type,rtcm->len,zcnt,staid,seqno);
    }
    if (type==3||type==22||type==23||type==24) {
        if (rtcm->staid!=0&&staid!=rtcm->staid) {
           trace(2,"rtcm2 station id changed: %d->%d\n",rtcm->staid,staid);
        }
        rtcm->staid=staid;
    }
    if (rtcm->staid!=0&&staid!=rtcm->staid) {
        trace(2,"rtcm2 station id invalid: %d %d\n",staid,rtcm->staid);
        return -1;
    }
    switch (type) {
        case  1: ret=od_rtk_rtcm2_decode_type1 (rtcm); break;
        case  3: ret=od_rtk_rtcm2_decode_type3 (rtcm); break;
        case  9: ret=od_rtk_rtcm2_decode_type1 (rtcm); break;
        case 14: ret=od_rtk_rtcm2_decode_type14(rtcm); break;
        case 16: ret=od_rtk_rtcm2_decode_type16(rtcm); break;
        case 17: ret=od_rtk_rtcm2_decode_type17(rtcm); break;
        case 18: ret=od_rtk_rtcm2_decode_type18(rtcm); break;
        case 19: ret=od_rtk_rtcm2_decode_type19(rtcm); break;
        case 22: ret=od_rtk_rtcm2_decode_type22(rtcm); break;
        case 23: ret=od_rtk_rtcm2_decode_type23(rtcm); break; /* not supported */
        case 24: ret=od_rtk_rtcm2_decode_type24(rtcm); break; /* not supported */
        case 31: ret=od_rtk_rtcm2_decode_type31(rtcm); break; /* not supported */
        case 32: ret=od_rtk_rtcm2_decode_type32(rtcm); break; /* not supported */
        case 34: ret=od_rtk_rtcm2_decode_type34(rtcm); break; /* not supported */
        case 36: ret=od_rtk_rtcm2_decode_type36(rtcm); break; /* not supported */
        case 37: ret=od_rtk_rtcm2_decode_type37(rtcm); break; /* not supported */
        case 59: ret=od_rtk_rtcm2_decode_type59(rtcm); break; /* not supported */
    }
    if (ret>=0) {
        if (1<=type&&type<=99) rtcm->nmsg2[type]++; else rtcm->nmsg2[0]++;
    }
    return ret;
}




/* ===== Embedded rtcm3.c ===== */
#pragma push_macro("P2_10")
#undef P2_10
#pragma push_macro("P2_28")
#undef P2_28
#pragma push_macro("P2_34")
#undef P2_34
#pragma push_macro("P2_41")
#undef P2_41
#pragma push_macro("P2_46")
#undef P2_46
#pragma push_macro("P2_59")
#undef P2_59
#pragma push_macro("P2_66")
#undef P2_66
#pragma push_macro("PRUNIT_GLO")
#undef PRUNIT_GLO
#pragma push_macro("PRUNIT_GPS")
#undef PRUNIT_GPS
#pragma push_macro("RANGE_MS")
#undef RANGE_MS

/*------------------------------------------------------------------------------
* rtcm3.c : RTCM ver.3 message decorder functions
*
*          Copyright (C) 2009-2020 by T.TAKASU, All rights reserved.
*
* references :
*     see rtcm.c
*
* version : $Revision:$ $Date:$
* history : 2012/05/14 1.0  separated from rtcm.c
*           2012/12/12 1.1  support gal/qzs ephemeris, gal/qzs ssr, msm
*                           add station id consistency test for obs data
*           2012/12/25 1.2  change compass msm id table
*           2013/01/31 1.3  change signal id by the latest draft (ref [13])
*           2013/02/23 1.4  change reference for rtcm 3 message (ref [14])
*           2013/05/19 1.5  gpst -> bdt of time-tag in beidou msm message
*           2014/05/02 1.6  fix bug on dropping last field of ssr message
*                           comply with rtcm 3.2 with amendment 1/2 (ref[15])
*                           delete MT 1046 according to ref [15]
*           2014/09/14 1.7  add receiver option -RT_INP
*           2014/12/06 1.8  support SBAS/BeiDou SSR messages (ref [16])
*           2015/03/22 1.9  add handling of iodcrc for beidou/sbas ssr messages
*           2015/04/27 1.10 support phase bias messages (MT2065-2070)
*           2015/09/07 1.11 add message count of MT 2000-2099
*           2015/10/21 1.12 add MT1046 support for IGS MGEX
*                           fix bug on decode of SSR 3/7 (code/phase bias)
*           2015/12/04 1.13 add MT63 beidou ephemeris (rtcm draft)
*                           fix bug on ssr 3 message decoding (#321)
*           2016/01/22 1.14 fix bug on L2C code in MT1004 (#131)
*           2016/08/20 1.15 fix bug on loss-of-lock detection in MSM 6/7 (#134)
*           2016/09/20 1.16 fix bug on MT1045 Galileo week rollover
*           2016/10/09 1.17 support MT1029 unicode text string
*           2017/04/11 1.18 fix bug on unchange-test of beidou ephemeris
*                           fix bug on week number in galileo ephemeris struct
*           2018/10/10 1.19 merge changes for 2.4.2 p13
*                           fix problem on eph.code for galileo ephemeris
*                           change mt for ssr 7 phase biases
*                           add rtcm option -GALINAV, -GALFNAV
*           2018/11/05 1.20 fix problem on invalid time in message monitor
*           2019/05/10 1.21 save galileo E5b data to obs index 2
*           2020/11/30 1.22 support MT1230 GLONASS code-phase biases
*                           support MT1131-1137,1041 (NavIC MSM and ephemeris)
*                           support MT4076 IGS SSR
*                           update MSM signal ID table (ref [17])
*                           update SSR signal and tracking mode ID table
*                           add week adjustment in MT1019,1044,1045,1046,1042
*                           use API code2idx() to get freq-index
*                           use API code2freq() to get carrier frequency
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/

/* constants -----------------------------------------------------------------*/

#define PRUNIT_GPS  299792.458  /* rtcm ver.3 unit of gps pseudorange (m) */
#define PRUNIT_GLO  599584.916  /* rtcm ver.3 unit of glonass pseudorange (m) */
#define RANGE_MS    (CLIGHT*0.001)      /* range in 1 ms */

#define P2_10       0.0009765625          /* 2^-10 */
#define P2_28       3.725290298461914E-09 /* 2^-28 */
#define P2_34       5.820766091346740E-11 /* 2^-34 */
#define P2_41       4.547473508864641E-13 /* 2^-41 */
#define P2_46       1.421085471520200E-14 /* 2^-46 */
#define P2_59       1.734723475976810E-18 /* 2^-59 */
#define P2_66       1.355252715606880E-20 /* 2^-66 */

/* type definition -----------------------------------------------------------*/

typedef struct {              /* multi-signal-message header type */
    uint8_t iod;              /* issue of data station */
    uint8_t time_s;           /* cumulative session transmitting time */
    uint8_t clk_str;          /* clock steering indicator */
    uint8_t clk_ext;          /* external clock indicator */
    uint8_t smooth;           /* divergence free smoothing indicator */
    uint8_t tint_s;           /* soothing interval */
    uint8_t nsat,nsig;        /* number of satellites/signals */
    uint8_t sats[64];         /* satellites */
    uint8_t sigs[32];         /* signals */
    uint8_t cellmask[64];     /* cell mask */
} od_rtk_rtcm3_msm_h_t;

/* MSM signal ID table -------------------------------------------------------*/
const char *msm_sig_gps[32]={
    /* GPS: ref [17] table 3.5-91 */
    ""  ,"1C","1P","1W",""  ,""  ,""  ,"2C","2P","2W",""  ,""  , /*  1-12 */
    ""  ,""  ,"2S","2L","2X",""  ,""  ,""  ,""  ,"5I","5Q","5X", /* 13-24 */
    ""  ,""  ,""  ,""  ,""  ,"1S","1L","1X"                      /* 25-32 */
};
const char *msm_sig_glo[32]={
    /* GLONASS: ref [17] table 3.5-96 */
    ""  ,"1C","1P",""  ,""  ,""  ,""  ,"2C","2P",""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_gal[32]={
    /* Galileo: ref [17] table 3.5-99 */
    ""  ,"1C","1A","1B","1X","1Z",""  ,"6C","6A","6B","6X","6Z",
    ""  ,"7I","7Q","7X",""  ,"8I","8Q","8X",""  ,"5I","5Q","5X",
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_qzs[32]={
    /* QZSS: ref [17] table 3.5-105 */
    ""  ,"1C",""  ,""  ,""  ,""  ,""  ,""  ,"6S","6L","6X",""  ,
    ""  ,""  ,"2S","2L","2X",""  ,""  ,""  ,""  ,"5I","5Q","5X",
    ""  ,""  ,""  ,""  ,""  ,"1S","1L","1X"
};
const char *msm_sig_sbs[32]={
    /* SBAS: ref [17] table 3.5-102 */
    ""  ,"1C",""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,"5I","5Q","5X",
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_cmp[32]={
    /* BeiDou: ref [17] table 3.5-108 */
    ""  ,"2I","2Q","2X",""  ,""  ,""  ,"6I","6Q","6X",""  ,""  ,
    ""  ,"7I","7Q","7X",""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
const char *msm_sig_irn[32]={
    /* NavIC/IRNSS: ref [17] table 3.5-108.3 */
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,""  ,"5A",""  ,""  ,
    ""  ,""  ,""  ,""  ,""  ,""  ,""  ,""
};
/* SSR signal and tracking mode IDs ------------------------------------------*/
const uint8_t ssr_sig_gps[32]={
    CODE_L1C,CODE_L1P,CODE_L1W,CODE_L1S,CODE_L1L,CODE_L2C,CODE_L2D,CODE_L2S,
    CODE_L2L,CODE_L2X,CODE_L2P,CODE_L2W,       0,       0,CODE_L5I,CODE_L5Q
};
const uint8_t ssr_sig_glo[32]={
    CODE_L1C,CODE_L1P,CODE_L2C,CODE_L2P,CODE_L4A,CODE_L4B,CODE_L6A,CODE_L6B,
    CODE_L3I,CODE_L3Q
};
const uint8_t ssr_sig_gal[32]={
    CODE_L1A,CODE_L1B,CODE_L1C,       0,       0,CODE_L5I,CODE_L5Q,       0,
    CODE_L7I,CODE_L7Q,       0,CODE_L8I,CODE_L8Q,       0,CODE_L6A,CODE_L6B,
    CODE_L6C
};
const uint8_t ssr_sig_qzs[32]={
    CODE_L1C,CODE_L1S,CODE_L1L,CODE_L2S,CODE_L2L,       0,CODE_L5I,CODE_L5Q,
           0,CODE_L6S,CODE_L6L,       0,       0,       0,       0,       0,
           0,CODE_L6E
};
const uint8_t ssr_sig_cmp[32]={
    CODE_L2I,CODE_L2Q,       0,CODE_L6I,CODE_L6Q,       0,CODE_L7I,CODE_L7Q,
           0,CODE_L1D,CODE_L1P,       0,CODE_L5D,CODE_L5P,       0,CODE_L1A,
           0,       0,CODE_L6A
};
const uint8_t ssr_sig_sbs[32]={
    CODE_L1C,CODE_L5I,CODE_L5Q
};
/* SSR update intervals ------------------------------------------------------*/
static const double od_rtk_rtcm3_ssrudint[16]={
    1,2,5,10,15,30,60,120,240,300,600,900,1800,3600,7200,10800
};
/* get sign-magnitude bits ---------------------------------------------------*/
static double od_rtk_rtcm3_getbitg(const uint8_t *buff, int pos, int len)
{
    double value=getbitu(buff,pos+1,len-1);
    return getbitu(buff,pos,1)?-value:value;
}
/* adjust weekly rollover of GPS time ----------------------------------------*/
static void od_rtk_rtcm3_adjweek(rtcm_t *rtcm, double tow)
{
    double tow_p;
    int week;
    
    /* if no time, get cpu time */
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tow_p=time2gpst(rtcm->time,&week);
    if      (tow<tow_p-302400.0) tow+=604800.0;
    else if (tow>tow_p+302400.0) tow-=604800.0;
    rtcm->time=gpst2time(week,tow);
}
/* adjust weekly rollover of BDS time ----------------------------------------*/
static int od_rtk_rtcm3_adjbdtweek(int week)
{
    int w;
    (void)time2bdt(gpst2bdt(utc2gpst(timeget())),&w);
    if (w<1) w=1; /* use 2006/1/1 if time is earlier than 2006/1/1 */
    return week+(w-week+512)/1024*1024;
}
/* adjust daily rollover of GLONASS time -------------------------------------*/
static void od_rtk_rtcm3_adjday_glot(rtcm_t *rtcm, double tod)
{
    gtime_t time;
    double tow,tod_p;
    int week;
    
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    time=timeadd(gpst2utc(rtcm->time),10800.0); /* glonass time */
    tow=time2gpst(time,&week);
    tod_p=fmod(tow,86400.0); tow-=tod_p;
    if      (tod<tod_p-43200.0) tod+=86400.0;
    else if (tod>tod_p+43200.0) tod-=86400.0;
    time=gpst2time(week,tow+tod);
    rtcm->time=utc2gpst(timeadd(time,-10800.0));
}
/* adjust carrier-phase rollover ---------------------------------------------*/
static double od_rtk_rtcm3_adjcp(rtcm_t *rtcm, int sat, int idx, double cp)
{
    if (rtcm->cp[sat-1][idx]==0.0) ;
    else if (cp<rtcm->cp[sat-1][idx]-750.0) cp+=1500.0;
    else if (cp>rtcm->cp[sat-1][idx]+750.0) cp-=1500.0;
    rtcm->cp[sat-1][idx]=cp;
    return cp;
}
/* loss-of-lock indicator ----------------------------------------------------*/
static int od_rtk_rtcm3_lossoflock(rtcm_t *rtcm, int sat, int idx, int lock)
{
    int lli=(!lock&&!rtcm->lock[sat-1][idx])||lock<rtcm->lock[sat-1][idx];
    rtcm->lock[sat-1][idx]=(uint16_t)lock;
    return lli;
}
/* S/N ratio -----------------------------------------------------------------*/
static uint16_t od_rtk_rtcm3_snratio(double snr)
{
    return (uint16_t)(snr<=0.0||100.0<=snr?0.0:snr/SNR_UNIT+0.5);
}
/* get observation data index ------------------------------------------------*/
static int od_rtk_rtcm3_obsindex(obs_t *obs, gtime_t time, int sat)
{
    int i,j;
    
    for (i=0;i<obs->n;i++) {
        if (obs->data[i].sat==sat) return i; /* field already exists */
    }
    if (i>=MAXOBS) return -1; /* overflow */
    
    /* add new field */
    obs->data[i].time=time;
    obs->data[i].sat=sat;
    for (j=0;j<NFREQ+NEXOBS;j++) {
        obs->data[i].L[j]=obs->data[i].P[j]=0.0;
        obs->data[i].D[j]=0.0;
        obs->data[i].SNR[j]=obs->data[i].LLI[j]=obs->data[i].code[j]=0;
    }
    obs->n++;
    return i;
}
/* test station ID consistency -----------------------------------------------*/
static int od_rtk_rtcm3_test_staid(rtcm_t *rtcm, int staid)
{
    char *p;
    int type,id;
    
    /* test station id option */
    if ((p=strstr(rtcm->opt,"-STA="))&&sscanf(p,"-STA=%d",&id)==1) {
        if (staid!=id) return 0;
    }
    /* save station id */
    if (rtcm->staid==0||rtcm->obsflag) {
        rtcm->staid=staid;
    }
    else if (staid!=rtcm->staid) {
        type=getbitu(rtcm->buff,24,12);
        trace(2,"rtcm3 %d staid invalid id=%d %d\n",type,staid,rtcm->staid);
        
        /* reset station id if station id error */
        rtcm->staid=0;
        return 0;
    }
    return 1;
}
/* decode type 1001-1004 message header --------------------------------------*/
static int od_rtk_rtcm3_decode_head1001(rtcm_t *rtcm, int *sync)
{
    double tow;
    char *msg,tstr[64];
    int i=24,staid,nsat,type;
    
    type=getbitu(rtcm->buff,i,12); i+=12;
    
    if (i+52<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12);       i+=12;
        tow  =getbitu(rtcm->buff,i,30)*0.001; i+=30;
        *sync=getbitu(rtcm->buff,i, 1);       i+= 1;
        nsat =getbitu(rtcm->buff,i, 5);
    }
    else {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    /* test station ID */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    od_rtk_rtcm3_adjweek(rtcm,tow);
    
    time2str(rtcm->time,tstr,2);
    trace(4,"decode_head1001: time=%s nsat=%d sync=%d\n",tstr,nsat,*sync);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d %s nsat=%2d sync=%d",staid,tstr,nsat,*sync);
    }
    return nsat;
}
/* decode type 1001: L1-only GPS RTK observation -----------------------------*/
static int od_rtk_rtcm3_decode_type1001(rtcm_t *rtcm)
{
    int sync;
    if (od_rtk_rtcm3_decode_head1001(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 1002: extended L1-only GPS RTK observables --------------------*/
static int od_rtk_rtcm3_decode_type1002(rtcm_t *rtcm)
{
    double pr1,cnr1,tt,cp1,freq=FREQ1;
    int i=24+64,j,index,nsat,sync,prn,code,sat,ppr1,lock1,amb,sys;
    
    if ((nsat=od_rtk_rtcm3_decode_head1001(rtcm,&sync))<0) return -1;
    
    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+74<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code =getbitu(rtcm->buff,i, 1); i+= 1;
        pr1  =getbitu(rtcm->buff,i,24); i+=24;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 8); i+= 8;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (prn<40) {
            sys=SYS_GPS;
        }
        else {
            sys=SYS_SBS; prn+=80;
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1002 satellite number error: prn=%d\n",prn);
            continue;
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=od_rtk_rtcm3_obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GPS;
        rtcm->obs.data[index].P[0]=pr1;
        
        if (ppr1!=(int)0xFFF80000) {
            cp1=od_rtk_rtcm3_adjcp(rtcm,sat,0,ppr1*0.0005*freq/CLIGHT);
            rtcm->obs.data[index].L[0]=pr1*freq/CLIGHT+cp1;
        }
        rtcm->obs.data[index].LLI[0]=od_rtk_rtcm3_lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=od_rtk_rtcm3_snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code?CODE_L1P:CODE_L1C;
    }
    return sync?0:1;
}
/* decode type 1003: L1&L2 gps rtk observables -------------------------------*/
static int od_rtk_rtcm3_decode_type1003(rtcm_t *rtcm)
{
    int sync;
    if (od_rtk_rtcm3_decode_head1001(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 1004: extended L1&L2 GPS RTK observables ----------------------*/
static int od_rtk_rtcm3_decode_type1004(rtcm_t *rtcm)
{
    const int L2codes[]={CODE_L2X,CODE_L2P,CODE_L2D,CODE_L2W};
    double pr1,cnr1,cnr2,tt,cp1,cp2,freq[2]={FREQ1,FREQ2};
    int i=24+64,j,index,nsat,sync,prn,sat,code1,code2,pr21,ppr1,ppr2;
    int lock1,lock2,amb,sys;
    
    if ((nsat=od_rtk_rtcm3_decode_head1001(rtcm,&sync))<0) return -1;
    
    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+125<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code1=getbitu(rtcm->buff,i, 1); i+= 1;
        pr1  =getbitu(rtcm->buff,i,24); i+=24;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 8); i+= 8;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        code2=getbitu(rtcm->buff,i, 2); i+= 2;
        pr21 =getbits(rtcm->buff,i,14); i+=14;
        ppr2 =getbits(rtcm->buff,i,20); i+=20;
        lock2=getbitu(rtcm->buff,i, 7); i+= 7;
        cnr2 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (prn<40) {
            sys=SYS_GPS;
        }
        else {
            sys=SYS_SBS; prn+=80;
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1004 satellite number error: sys=%d prn=%d\n",sys,prn);
            continue;
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=od_rtk_rtcm3_obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GPS;
        rtcm->obs.data[index].P[0]=pr1;
        
        if (ppr1!=(int)0xFFF80000) {
            cp1=od_rtk_rtcm3_adjcp(rtcm,sat,0,ppr1*0.0005*freq[0]/CLIGHT);
            rtcm->obs.data[index].L[0]=pr1*freq[0]/CLIGHT+cp1;
        }
        rtcm->obs.data[index].LLI[0]=od_rtk_rtcm3_lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=od_rtk_rtcm3_snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code1?CODE_L1P:CODE_L1C;
        
        if (pr21!=(int)0xFFFFE000) {
            rtcm->obs.data[index].P[1]=pr1+pr21*0.02;
        }
        if (ppr2!=(int)0xFFF80000) {
            cp2=od_rtk_rtcm3_adjcp(rtcm,sat,1,ppr2*0.0005*freq[1]/CLIGHT);
            rtcm->obs.data[index].L[1]=pr1*freq[1]/CLIGHT+cp2;
        }
        rtcm->obs.data[index].LLI[1]=od_rtk_rtcm3_lossoflock(rtcm,sat,1,lock2);
        rtcm->obs.data[index].SNR[1]=od_rtk_rtcm3_snratio(cnr2*0.25);
        rtcm->obs.data[index].code[1]=L2codes[code2];
    }
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* get signed 38bit field ----------------------------------------------------*/
static double od_rtk_rtcm3_getbits_38(const uint8_t *buff, int pos)
{
    return (double)getbits(buff,pos,32)*64.0+getbitu(buff,pos+32,6);
}
/* decode type 1005: stationary RTK reference station ARP --------------------*/
static int od_rtk_rtcm3_decode_type1005(rtcm_t *rtcm)
{
    double rr[3],re[3],pos[3];
    char *msg;
    int i=24+12,j,staid,itrf;
    
    if (i+140==rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12;
        itrf =getbitu(rtcm->buff,i, 6); i+= 6+4;
        rr[0]=od_rtk_rtcm3_getbits_38(rtcm->buff,i); i+=38+2;
        rr[1]=od_rtk_rtcm3_getbits_38(rtcm->buff,i); i+=38+2;
        rr[2]=od_rtk_rtcm3_getbits_38(rtcm->buff,i);
    }
    else {
        trace(2,"rtcm3 1005 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        for (j=0;j<3;j++) re[j]=rr[j]*0.0001;
        ecef2pos(re,pos);
        sprintf(msg," staid=%4d pos=%.8f %.8f %.3f",staid,pos[0]*R2D,pos[1]*R2D,
                pos[2]);
    }
    /* test station id */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    sprintf(rtcm->sta.name,"%04d",staid);
    rtcm->sta.deltype=0; /* xyz */
    for (j=0;j<3;j++) {
        rtcm->sta.pos[j]=rr[j]*0.0001;
        rtcm->sta.del[j]=0.0;
    }
    rtcm->sta.hgt=0.0;
    rtcm->sta.itrf=itrf;
    return 5;
}
/* decode type 1006: stationary RTK reference station ARP with height --------*/
static int od_rtk_rtcm3_decode_type1006(rtcm_t *rtcm)
{
    double rr[3],re[3],pos[3],anth;
    char *msg;
    int i=24+12,j,staid,itrf;
    
    if (i+156<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12;
        itrf =getbitu(rtcm->buff,i, 6); i+= 6+4;
        rr[0]=od_rtk_rtcm3_getbits_38(rtcm->buff,i); i+=38+2;
        rr[1]=od_rtk_rtcm3_getbits_38(rtcm->buff,i); i+=38+2;
        rr[2]=od_rtk_rtcm3_getbits_38(rtcm->buff,i); i+=38;
        anth =getbitu(rtcm->buff,i,16);
    }
    else {
        trace(2,"rtcm3 1006 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        for (j=0;j<3;j++) re[j]=rr[j]*0.0001;
        ecef2pos(re,pos);
        sprintf(msg," staid=%4d pos=%.8f %.8f %.3f anth=%.3f",staid,pos[0]*R2D,
                pos[1]*R2D,pos[2],anth*0.0001);
    }
    /* test station id */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    sprintf(rtcm->sta.name,"%04d",staid);
    rtcm->sta.deltype=1; /* xyz */
    for (j=0;j<3;j++) {
        rtcm->sta.pos[j]=rr[j]*0.0001;
        rtcm->sta.del[j]=0.0;
    }
    rtcm->sta.hgt=anth*0.0001;
    rtcm->sta.itrf=itrf;
    return 5;
}
/* decode type 1007: antenna descriptor --------------------------------------*/
static int od_rtk_rtcm3_decode_type1007(rtcm_t *rtcm)
{
    char des[32]="";
    char *msg;
    int i=24+12,j,staid,n,setup;
    
    n=getbitu(rtcm->buff,i+12,8);
    
    if (i+28+8*n<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12+8;
        for (j=0;j<n&&j<31;j++) {
            des[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        setup=getbitu(rtcm->buff,i, 8);
    }
    else {
        trace(2,"rtcm3 1007 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d",staid);
    }
    /* test station ID */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    sprintf(rtcm->sta.name,"%04d",staid);
    strncpy(rtcm->sta.antdes,des,n); rtcm->sta.antdes[n]='\0';
    rtcm->sta.antsetup=setup;
    rtcm->sta.antsno[0]='\0';
    return 5;
}
/* decode type 1008: antenna descriptor & serial number ----------------------*/
static int od_rtk_rtcm3_decode_type1008(rtcm_t *rtcm)
{
    char des[32]="",sno[32]="";
    char *msg;
    int i=24+12,j,staid,n,m,setup;
    
    n=getbitu(rtcm->buff,i+12,8);
    m=getbitu(rtcm->buff,i+28+8*n,8);
    
    if (i+36+8*(n+m)<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12+8;
        for (j=0;j<n&&j<31;j++) {
            des[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        setup=getbitu(rtcm->buff,i, 8); i+=8+8;
        for (j=0;j<m&&j<31;j++) {
            sno[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
    }
    else {
        trace(2,"rtcm3 1008 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d",staid);
    }
    /* test station ID */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    sprintf(rtcm->sta.name,"%04d",staid);
    strncpy(rtcm->sta.antdes,des,n); rtcm->sta.antdes[n]='\0';
    rtcm->sta.antsetup=setup;
    strncpy(rtcm->sta.antsno,sno,m); rtcm->sta.antsno[m]='\0';
    return 5;
}
/* decode type 1009-1012 message header --------------------------------------*/
static int od_rtk_rtcm3_decode_head1009(rtcm_t *rtcm, int *sync)
{
    double tod;
    char *msg,tstr[64];
    int i=24,staid,nsat,type;
    
    type=getbitu(rtcm->buff,i,12); i+=12;
    
    if (i+49<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12);       i+=12;
        tod  =getbitu(rtcm->buff,i,27)*0.001; i+=27; /* sec in a day */
        *sync=getbitu(rtcm->buff,i, 1);       i+= 1;
        nsat =getbitu(rtcm->buff,i, 5);
    }
    else {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    /* test station ID */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    od_rtk_rtcm3_adjday_glot(rtcm,tod);
    
    time2str(rtcm->time,tstr,2);
    trace(4,"decode_head1009: time=%s nsat=%d sync=%d\n",tstr,nsat,*sync);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d %s nsat=%2d sync=%d",staid,tstr,nsat,*sync);
    }
    return nsat;
}
/* decode type 1009: L1-only glonass rtk observables -------------------------*/
static int od_rtk_rtcm3_decode_type1009(rtcm_t *rtcm)
{
    int sync;
    if (od_rtk_rtcm3_decode_head1009(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 1010: extended L1-only glonass rtk observables ----------------*/
static int od_rtk_rtcm3_decode_type1010(rtcm_t *rtcm)
{
    double pr1,cnr1,tt,cp1,freq1;
    int i=24+61,j,index,nsat,sync,prn,sat,code,fcn,ppr1,lock1,amb,sys=SYS_GLO;
    
    if ((nsat=od_rtk_rtcm3_decode_head1009(rtcm,&sync))<0) return -1;
    
    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+79<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code =getbitu(rtcm->buff,i, 1); i+= 1;
        fcn  =getbitu(rtcm->buff,i, 5); i+= 5; /* fcn+7 */
        pr1  =getbitu(rtcm->buff,i,25); i+=25;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 7); i+= 7;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1010 satellite number error: prn=%d\n",prn);
            continue;
        }
        if (!rtcm->nav.glo_fcn[prn-1]) {
            rtcm->nav.glo_fcn[prn-1]=fcn-7+8; /* fcn+8 */
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=od_rtk_rtcm3_obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GLO;
        rtcm->obs.data[index].P[0]=pr1;
        
        if (ppr1!=(int)0xFFF80000) {
            freq1=code2freq(SYS_GLO,CODE_L1C,fcn-7);
            cp1=od_rtk_rtcm3_adjcp(rtcm,sat,0,ppr1*0.0005*freq1/CLIGHT);
            rtcm->obs.data[index].L[0]=pr1*freq1/CLIGHT+cp1;
        }
        rtcm->obs.data[index].LLI[0]=od_rtk_rtcm3_lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=od_rtk_rtcm3_snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code?CODE_L1P:CODE_L1C;
    }
    return sync?0:1;
}
/* decode type 1011: L1&L2 GLONASS RTK observables ---------------------------*/
static int od_rtk_rtcm3_decode_type1011(rtcm_t *rtcm)
{
    int sync;
    if (od_rtk_rtcm3_decode_head1009(rtcm,&sync)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 1012: extended L1&L2 GLONASS RTK observables ------------------*/
static int od_rtk_rtcm3_decode_type1012(rtcm_t *rtcm)
{
    double pr1,cnr1,cnr2,tt,cp1,cp2,freq1,freq2;
    int i=24+61,j,index,nsat,sync,prn,sat,fcn,code1,code2,pr21,ppr1,ppr2;
    int lock1,lock2,amb,sys=SYS_GLO;
    
    if ((nsat=od_rtk_rtcm3_decode_head1009(rtcm,&sync))<0) return -1;
    
    for (j=0;j<nsat&&rtcm->obs.n<MAXOBS&&i+130<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i, 6); i+= 6;
        code1=getbitu(rtcm->buff,i, 1); i+= 1;
        fcn  =getbitu(rtcm->buff,i, 5); i+= 5; /* fcn+7 */
        pr1  =getbitu(rtcm->buff,i,25); i+=25;
        ppr1 =getbits(rtcm->buff,i,20); i+=20;
        lock1=getbitu(rtcm->buff,i, 7); i+= 7;
        amb  =getbitu(rtcm->buff,i, 7); i+= 7;
        cnr1 =getbitu(rtcm->buff,i, 8); i+= 8;
        code2=getbitu(rtcm->buff,i, 2); i+= 2;
        pr21 =getbits(rtcm->buff,i,14); i+=14;
        ppr2 =getbits(rtcm->buff,i,20); i+=20;
        lock2=getbitu(rtcm->buff,i, 7); i+= 7;
        cnr2 =getbitu(rtcm->buff,i, 8); i+= 8;
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 1012 satellite number error: sys=%d prn=%d\n",sys,prn);
            continue;
        }
        if (!rtcm->nav.glo_fcn[prn-1]) {
            rtcm->nav.glo_fcn[prn-1]=fcn-7+8; /* fcn+8 */
        }
        tt=timediff(rtcm->obs.data[0].time,rtcm->time);
        if (rtcm->obsflag||fabs(tt)>1E-9) {
            rtcm->obs.n=rtcm->obsflag=0;
        }
        if ((index=od_rtk_rtcm3_obsindex(&rtcm->obs,rtcm->time,sat))<0) continue;
        pr1=pr1*0.02+amb*PRUNIT_GLO;
        rtcm->obs.data[index].P[0]=pr1;
        
        if (ppr1!=(int)0xFFF80000) {
            freq1=code2freq(SYS_GLO,CODE_L1C,fcn-7);
            cp1=od_rtk_rtcm3_adjcp(rtcm,sat,0,ppr1*0.0005*freq1/CLIGHT);
            rtcm->obs.data[index].L[0]=pr1*freq1/CLIGHT+cp1;
        }
        rtcm->obs.data[index].LLI[0]=od_rtk_rtcm3_lossoflock(rtcm,sat,0,lock1);
        rtcm->obs.data[index].SNR[0]=od_rtk_rtcm3_snratio(cnr1*0.25);
        rtcm->obs.data[index].code[0]=code1?CODE_L1P:CODE_L1C;
        
        if (pr21!=(int)0xFFFFE000) {
            rtcm->obs.data[index].P[1]=pr1+pr21*0.02;
        }
        if (ppr2!=(int)0xFFF80000) {
            freq2=code2freq(SYS_GLO,CODE_L2C,fcn-7);
            cp2=od_rtk_rtcm3_adjcp(rtcm,sat,1,ppr2*0.0005*freq2/CLIGHT);
            rtcm->obs.data[index].L[1]=pr1*freq2/CLIGHT+cp2;
        }
        rtcm->obs.data[index].LLI[1]=od_rtk_rtcm3_lossoflock(rtcm,sat,1,lock2);
        rtcm->obs.data[index].SNR[1]=od_rtk_rtcm3_snratio(cnr2*0.25);
        rtcm->obs.data[index].code[1]=code2?CODE_L2P:CODE_L2C;
    }
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 1013: system parameters ---------------------------------------*/
static int od_rtk_rtcm3_decode_type1013(rtcm_t *rtcm)
{
    return 0;
}
/* decode type 1019: GPS ephemerides -----------------------------------------*/
static int od_rtk_rtcm3_decode_type1019(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA,tt;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_GPS;
    
    if (i+476<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.code  =getbitu(rtcm->buff,i, 2);              i+= 2;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        eph.iodc  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.svh   =getbitu(rtcm->buff,i, 6);              i+= 6;
        eph.flag  =getbitu(rtcm->buff,i, 1);              i+= 1;
        eph.fit   =getbitu(rtcm->buff,i, 1)?0.0:4.0; /* 0:4hr,1:>4hr */
    }
    else {
        trace(2,"rtcm3 1019 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (prn>=40) {
        sys=SYS_SBS; prn+=80;
    }
    trace(4,"decode_type1019: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
                prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1019 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjgpsweek(week);
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tt=timediff(gpst2time(eph.week,eph.toes),rtcm->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt>=302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[sat-1].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[sat-1]=eph;
    rtcm->ephsat=sat;
    rtcm->ephset=0;
    return 2;
}
/* decode type 1020: GLONASS ephemerides -------------------------------------*/
static int od_rtk_rtcm3_decode_type1020(rtcm_t *rtcm)
{
    geph_t geph={0};
    double tk_h,tk_m,tk_s,toe,tow,tod,tof;
    char *msg;
    int i=24+12,prn,sat,week,tb,bn,sys=SYS_GLO;
    
    if (i+348<=rtcm->len*8) {
        prn        =getbitu(rtcm->buff,i, 6);           i+= 6;
        geph.frq   =getbitu(rtcm->buff,i, 5)-7;         i+= 5+2+2;
        tk_h       =getbitu(rtcm->buff,i, 5);           i+= 5;
        tk_m       =getbitu(rtcm->buff,i, 6);           i+= 6;
        tk_s       =getbitu(rtcm->buff,i, 1)*30.0;      i+= 1;
        bn         =getbitu(rtcm->buff,i, 1);           i+= 1+1;
        tb         =getbitu(rtcm->buff,i, 7);           i+= 7;
        geph.vel[0]=od_rtk_rtcm3_getbitg(rtcm->buff,i,24)*P2_20*1E3; i+=24;
        geph.pos[0]=od_rtk_rtcm3_getbitg(rtcm->buff,i,27)*P2_11*1E3; i+=27;
        geph.acc[0]=od_rtk_rtcm3_getbitg(rtcm->buff,i, 5)*P2_30*1E3; i+= 5;
        geph.vel[1]=od_rtk_rtcm3_getbitg(rtcm->buff,i,24)*P2_20*1E3; i+=24;
        geph.pos[1]=od_rtk_rtcm3_getbitg(rtcm->buff,i,27)*P2_11*1E3; i+=27;
        geph.acc[1]=od_rtk_rtcm3_getbitg(rtcm->buff,i, 5)*P2_30*1E3; i+= 5;
        geph.vel[2]=od_rtk_rtcm3_getbitg(rtcm->buff,i,24)*P2_20*1E3; i+=24;
        geph.pos[2]=od_rtk_rtcm3_getbitg(rtcm->buff,i,27)*P2_11*1E3; i+=27;
        geph.acc[2]=od_rtk_rtcm3_getbitg(rtcm->buff,i, 5)*P2_30*1E3; i+= 5+1;
        geph.gamn  =od_rtk_rtcm3_getbitg(rtcm->buff,i,11)*P2_40;     i+=11+3;
        geph.taun  =od_rtk_rtcm3_getbitg(rtcm->buff,i,22)*P2_30;     i+=22;
        geph.dtaun =od_rtk_rtcm3_getbitg(rtcm->buff,i, 5)*P2_30;     i+=5;
        geph.age   =getbitu(rtcm->buff,i, 5);
    }
    else {
        trace(2,"rtcm3 1020 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1020 satellite number error: prn=%d\n",prn);
        return -1;
    }
    trace(4,"decode_type1020: prn=%d tk=%02.0f:%02.0f:%02.0f\n",prn,tk_h,tk_m,tk_s);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d tk=%02.0f:%02.0f:%02.0f frq=%2d bn=%d tb=%d",
                prn,tk_h,tk_m,tk_s,geph.frq,bn,tb);
    }
    geph.sat=sat;
    geph.svh=bn;
    geph.iode=tb&0x7F;
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tow=time2gpst(gpst2utc(rtcm->time),&week);
    tod=fmod(tow,86400.0); tow-=tod;
    tof=tk_h*3600.0+tk_m*60.0+tk_s-10800.0; /* lt->utc */
    if      (tof<tod-43200.0) tof+=86400.0;
    else if (tof>tod+43200.0) tof-=86400.0;
    geph.tof=utc2gpst(gpst2time(week,tow+tof));
    toe=tb*900.0-10800.0; /* lt->utc */
    if      (toe<tod-43200.0) toe+=86400.0;
    else if (toe>tod+43200.0) toe-=86400.0;
    geph.toe=utc2gpst(gpst2time(week,tow+toe)); /* utc->gpst */
    
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (fabs(timediff(geph.toe,rtcm->nav.geph[prn-1].toe))<1.0&&
            geph.svh==rtcm->nav.geph[prn-1].svh) return 0; /* unchanged */
    }
    rtcm->nav.geph[prn-1]=geph;
    rtcm->ephsat=sat;
    rtcm->ephset=0;
    return 2;
}
/* decode type 1021: helmert/abridged molodenski -----------------------------*/
static int od_rtk_rtcm3_decode_type1021(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1021: not supported message\n");
    return 0;
}
/* decode type 1022: Moledenski-Badekas transfromation -----------------------*/
static int od_rtk_rtcm3_decode_type1022(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1022: not supported message\n");
    return 0;
}
/* decode type 1023: residual, ellipsoidal grid representation ---------------*/
static int od_rtk_rtcm3_decode_type1023(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1023: not supported message\n");
    return 0;
}
/* decode type 1024: residual, plane grid representation ---------------------*/
static int od_rtk_rtcm3_decode_type1024(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1024: not supported message\n");
    return 0;
}
/* decode type 1025: projection (types except LCC2SP,OM) ---------------------*/
static int od_rtk_rtcm3_decode_type1025(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1025: not supported message\n");
    return 0;
}
/* decode type 1026: projection (LCC2SP - lambert conic conformal (2sp)) -----*/
static int od_rtk_rtcm3_decode_type1026(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1026: not supported message\n");
    return 0;
}
/* decode type 1027: projection (type OM - oblique mercator) -----------------*/
static int od_rtk_rtcm3_decode_type1027(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1027: not supported message\n");
    return 0;
}
/* decode type 1029: UNICODE text string -------------------------------------*/
static int od_rtk_rtcm3_decode_type1029(rtcm_t *rtcm)
{
    char *msg;
    int i=24+12,j,staid,mjd,tod,nchar,cunit;
    
    if (i+60<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12;
        mjd  =getbitu(rtcm->buff,i,16); i+=16;
        tod  =getbitu(rtcm->buff,i,17); i+=17;
        nchar=getbitu(rtcm->buff,i, 7); i+= 7;
        cunit=getbitu(rtcm->buff,i, 8); i+= 8;
    }
    else {
        trace(2,"rtcm3 1029 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (i+nchar*8>rtcm->len*8) {
        trace(2,"rtcm3 1029 length error: len=%d nchar=%d\n",rtcm->len,nchar);
        return -1;
    } 
    for (j=0;j<nchar&&j<126;j++) {
        rtcm->msg[j]=getbitu(rtcm->buff,i,8); i+=8;
    }
    rtcm->msg[j]='\0';
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d text=%s",staid,rtcm->msg);
    }
    return 0;
}
/* decode type 1030: network RTK residual ------------------------------------*/
static int od_rtk_rtcm3_decode_type1030(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1030: not supported message\n");
    return 0;
}
/* decode type 1031: GLONASS network RTK residual ----------------------------*/
static int od_rtk_rtcm3_decode_type1031(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1031: not supported message\n");
    return 0;
}
/* decode type 1032: physical reference station position information ---------*/
static int od_rtk_rtcm3_decode_type1032(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1032: not supported message\n");
    return 0;
}
/* decode type 1033: receiver and antenna descriptor -------------------------*/
static int od_rtk_rtcm3_decode_type1033(rtcm_t *rtcm)
{
    char des[32]="",sno[32]="",rec[32]="",ver[32]="",rsn[32]="";
    char *msg;
    int i=24+12,j,staid,n,m,n1,n2,n3,setup;
    
    n =getbitu(rtcm->buff,i+12,8);
    m =getbitu(rtcm->buff,i+28+8*n,8);
    n1=getbitu(rtcm->buff,i+36+8*(n+m),8);
    n2=getbitu(rtcm->buff,i+44+8*(n+m+n1),8);
    n3=getbitu(rtcm->buff,i+52+8*(n+m+n1+n2),8);
    
    if (i+60+8*(n+m+n1+n2+n3)<=rtcm->len*8) {
        staid=getbitu(rtcm->buff,i,12); i+=12+8;
        for (j=0;j<n&&j<31;j++) {
            des[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        setup=getbitu(rtcm->buff,i, 8); i+=8+8;
        for (j=0;j<m&&j<31;j++) {
            sno[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        i+=8;
        for (j=0;j<n1&&j<31;j++) {
            rec[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        i+=8;
        for (j=0;j<n2&&j<31;j++) {
            ver[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
        i+=8;
        for (j=0;j<n3&&j<31;j++) {
            rsn[j]=(char)getbitu(rtcm->buff,i,8); i+=8;
        }
    }
    else {
        trace(2,"rtcm3 1033 length error: len=%d\n",rtcm->len);
        return -1;
    }
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d",staid);
    }
    /* test station id */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    sprintf(rtcm->sta.name,"%04d",staid);
    strncpy(rtcm->sta.antdes, des,n ); rtcm->sta.antdes [n] ='\0';
    rtcm->sta.antsetup=setup;
    strncpy(rtcm->sta.antsno, sno,m ); rtcm->sta.antsno [m] ='\0';
    strncpy(rtcm->sta.rectype,rec,n1); rtcm->sta.rectype[n1]='\0';
    strncpy(rtcm->sta.recver, ver,n2); rtcm->sta.recver [n2]='\0';
    strncpy(rtcm->sta.recsno, rsn,n3); rtcm->sta.recsno [n3]='\0';
    
    trace(3,"rtcm3 1033: ant=%s:%s rec=%s:%s:%s\n",des,sno,rec,ver,rsn);
    return 5;
}
/* decode type 1034: GPS network FKP gradient --------------------------------*/
static int od_rtk_rtcm3_decode_type1034(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1034: not supported message\n");
    return 0;
}
/* decode type 1035: GLONASS network FKP gradient ----------------------------*/
static int od_rtk_rtcm3_decode_type1035(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1035: not supported message\n");
    return 0;
}
/* decode type 1037: GLONASS network RTK ionospheric correction difference ---*/
static int od_rtk_rtcm3_decode_type1037(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1037: not supported message\n");
    return 0;
}
/* decode type 1038: GLONASS network RTK geometic correction difference ------*/
static int od_rtk_rtcm3_decode_type1038(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1038: not supported message\n");
    return 0;
}
/* decode type 1039: GLONASS network RTK combined correction difference ------*/
static int od_rtk_rtcm3_decode_type1039(rtcm_t *rtcm)
{
    trace(2,"rtcm3 1039: not supported message\n");
    return 0;
}
/* decode type 1041: NavIC/IRNSS ephemerides ---------------------------------*/
static int od_rtk_rtcm3_decode_type1041(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA,tt;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_IRN;
    
    if (i+482-12<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.deln  =getbits(rtcm->buff,i,22)*P2_41*SC2RAD; i+=22;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8+10; /* IODEC */
        eph.svh   =getbitu(rtcm->buff,i, 2);              i+= 2; /* L5+Sflag */
        eph.cuc   =getbits(rtcm->buff,i,15)*P2_28;        i+=15;
        eph.cus   =getbits(rtcm->buff,i,15)*P2_28;        i+=15;
        eph.cic   =getbits(rtcm->buff,i,15)*P2_28;        i+=15;
        eph.cis   =getbits(rtcm->buff,i,15)*P2_28;        i+=15;
        eph.crc   =getbits(rtcm->buff,i,15)*0.0625;       i+=15;
        eph.crs   =getbits(rtcm->buff,i,15)*0.0625;       i+=15;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,22)*P2_41*SC2RAD; i+=22;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD;
    }
    else {
        trace(2,"rtcm3 1041 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1041: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
                prn,eph.iode,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1041 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjgpsweek(week);
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tt=timediff(gpst2time(eph.week,eph.toes),rtcm->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt>=302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    eph.iodc=eph.iode;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[sat-1].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[sat-1]=eph;
    rtcm->ephsat=sat;
    rtcm->ephset=0;
    return 2;
}
/* decode type 1044: QZSS ephemerides ----------------------------------------*/
static int od_rtk_rtcm3_decode_type1044(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA,tt;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_QZS;
    
    if (i+473<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 4)+192;          i+= 4;
        toc       =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.f2    =getbits(rtcm->buff,i, 8)*P2_55;        i+= 8;
        eph.f1    =getbits(rtcm->buff,i,16)*P2_43;        i+=16;
        eph.f0    =getbits(rtcm->buff,i,22)*P2_31;        i+=22;
        eph.iode  =getbitu(rtcm->buff,i, 8);              i+= 8;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,16)*16.0;         i+=16;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.code  =getbitu(rtcm->buff,i, 2);              i+= 2;
        week      =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.svh   =getbitu(rtcm->buff,i, 6);              i+= 6;
        eph.tgd[0]=getbits(rtcm->buff,i, 8)*P2_31;        i+= 8;
        eph.iodc  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.fit   =getbitu(rtcm->buff,i, 1)?0.0:2.0; /* 0:2hr,1:>2hr */
    }
    else {
        trace(2,"rtcm3 1044 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1044: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%3d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
                prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1044 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=adjgpsweek(week);
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tt=timediff(gpst2time(eph.week,eph.toes),rtcm->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt>=302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    eph.flag=1; /* fixed to 1 */
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[sat-1].iode&&
            eph.iodc==rtcm->nav.eph[sat-1].iodc) return 0; /* unchanged */
    }
    rtcm->nav.eph[sat-1]=eph;
    rtcm->ephsat=sat;
    rtcm->ephset=0;
    return 2;
}
/* decode type 1045: Galileo F/NAV satellite ephemerides ---------------------*/
static int od_rtk_rtcm3_decode_type1045(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA,tt;
    char *msg;
    int i=24+12,prn,sat,week,e5a_hs,e5a_dvs,rsv,sys=SYS_GAL;
    
    if (strstr(rtcm->opt,"-GALINAV")) return 0;

    if (i+484<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,12);              i+=12; /* gst-week */
        eph.iode  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 8);              i+= 8;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        toc       =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.f2    =getbits(rtcm->buff,i, 6)*P2_59;        i+= 6;
        eph.f1    =getbits(rtcm->buff,i,21)*P2_46;        i+=21;
        eph.f0    =getbits(rtcm->buff,i,31)*P2_34;        i+=31;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i,10)*P2_32;        i+=10; /* E5a/E1 */
        e5a_hs    =getbitu(rtcm->buff,i, 2);              i+= 2; /* OSHS */
        e5a_dvs   =getbitu(rtcm->buff,i, 1);              i+= 1; /* OSDVS */
        rsv       =getbitu(rtcm->buff,i, 7);
    }
    else {
        trace(2,"rtcm3 1045 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1045: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d week=%d toe=%6.0f toc=%6.0f hs=%d dvs=%d",
                prn,eph.iode,week,eph.toes,toc,e5a_hs,e5a_dvs);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1045 satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (strstr(rtcm->opt,"-GALINAV")) {
        return 0;
    }
    eph.sat=sat;
    eph.week=week+1024; /* gal-week = gst-week + 1024 */
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tt=timediff(gpst2time(eph.week,eph.toes),rtcm->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt>=302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    eph.svh=(e5a_hs<<4)+(e5a_dvs<<3);
    eph.code=(1<<1)+(1<<8); /* data source = F/NAV+E5a */
    eph.iodc=eph.iode;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[sat-1+MAXSAT].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[sat-1+MAXSAT]=eph;
    rtcm->ephsat=sat;
    rtcm->ephset=1; /* F/NAV */
    return 2;
}
/* decode type 1046: Galileo I/NAV satellite ephemerides ---------------------*/
static int od_rtk_rtcm3_decode_type1046(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA,tt;
    char *msg;
    int i=24+12,prn,sat,week,e5b_hs,e5b_dvs,e1_hs,e1_dvs,sys=SYS_GAL;
    
    if (strstr(rtcm->opt,"-GALFNAV")) return 0;

    if (i+492<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,12);              i+=12;
        eph.iode  =getbitu(rtcm->buff,i,10);              i+=10;
        eph.sva   =getbitu(rtcm->buff,i, 8);              i+= 8;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        toc       =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.f2    =getbits(rtcm->buff,i, 6)*P2_59;        i+= 6;
        eph.f1    =getbits(rtcm->buff,i,21)*P2_46;        i+=21;
        eph.f0    =getbits(rtcm->buff,i,31)*P2_34;        i+=31;
        eph.crs   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,14)*60.0;         i+=14;
        eph.cic   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,16)*P2_29;        i+=16;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,16)*P2_5;         i+=16;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i,10)*P2_32;        i+=10; /* E5a/E1 */
        eph.tgd[1]=getbits(rtcm->buff,i,10)*P2_32;        i+=10; /* E5b/E1 */
        e5b_hs    =getbitu(rtcm->buff,i, 2);              i+= 2; /* E5b OSHS */
        e5b_dvs   =getbitu(rtcm->buff,i, 1);              i+= 1; /* E5b OSDVS */
        e1_hs     =getbitu(rtcm->buff,i, 2);              i+= 2; /* E1 OSHS */
        e1_dvs    =getbitu(rtcm->buff,i, 1);              i+= 1; /* E1 OSDVS */
    }
    else {
        trace(2,"rtcm3 1046 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1046: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d week=%d toe=%6.0f toc=%6.0f hs=%d %d dvs=%d %d",
                prn,eph.iode,week,eph.toes,toc,e5b_hs,e1_hs,e5b_dvs,e1_dvs);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1046 satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (strstr(rtcm->opt,"-GALFNAV")) {
        return 0;
    }
    eph.sat=sat;
    eph.week=week+1024; /* gal-week = gst-week + 1024 */
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tt=timediff(gpst2time(eph.week,eph.toes),rtcm->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt>=302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=gpst2time(eph.week,toc);
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    eph.svh=(e5b_hs<<7)+(e5b_dvs<<6)+(e1_hs<<1)+(e1_dvs<<0);
    eph.code=(1<<0)+(1<<2)+(1<<9); /* data source = I/NAV+E1+E5b */
    eph.iodc=eph.iode;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (eph.iode==rtcm->nav.eph[sat-1].iode) return 0; /* unchanged */
    }
    rtcm->nav.eph[sat-1]=eph;
    rtcm->ephsat=sat;
    rtcm->ephset=0; /* I/NAV */
    return 2;
}
/* decode type 1042/63: Beidou ephemerides -----------------------------------*/
static int od_rtk_rtcm3_decode_type1042(rtcm_t *rtcm)
{
    eph_t eph={0};
    double toc,sqrtA,tt;
    char *msg;
    int i=24+12,prn,sat,week,sys=SYS_CMP;
    
    if (i+499<=rtcm->len*8) {
        prn       =getbitu(rtcm->buff,i, 6);              i+= 6;
        week      =getbitu(rtcm->buff,i,13);              i+=13;
        eph.sva   =getbitu(rtcm->buff,i, 4);              i+= 4;
        eph.idot  =getbits(rtcm->buff,i,14)*P2_43*SC2RAD; i+=14;
        eph.iode  =getbitu(rtcm->buff,i, 5);              i+= 5; /* AODE */
        toc       =getbitu(rtcm->buff,i,17)*8.0;          i+=17;
        eph.f2    =getbits(rtcm->buff,i,11)*P2_66;        i+=11;
        eph.f1    =getbits(rtcm->buff,i,22)*P2_50;        i+=22;
        eph.f0    =getbits(rtcm->buff,i,24)*P2_33;        i+=24;
        eph.iodc  =getbitu(rtcm->buff,i, 5);              i+= 5; /* AODC */
        eph.crs   =getbits(rtcm->buff,i,18)*P2_6;         i+=18;
        eph.deln  =getbits(rtcm->buff,i,16)*P2_43*SC2RAD; i+=16;
        eph.M0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cuc   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        eph.e     =getbitu(rtcm->buff,i,32)*P2_33;        i+=32;
        eph.cus   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        sqrtA     =getbitu(rtcm->buff,i,32)*P2_19;        i+=32;
        eph.toes  =getbitu(rtcm->buff,i,17)*8.0;          i+=17;
        eph.cic   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        eph.OMG0  =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.cis   =getbits(rtcm->buff,i,18)*P2_31;        i+=18;
        eph.i0    =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.crc   =getbits(rtcm->buff,i,18)*P2_6;         i+=18;
        eph.omg   =getbits(rtcm->buff,i,32)*P2_31*SC2RAD; i+=32;
        eph.OMGd  =getbits(rtcm->buff,i,24)*P2_43*SC2RAD; i+=24;
        eph.tgd[0]=getbits(rtcm->buff,i,10)*1E-10;        i+=10;
        eph.tgd[1]=getbits(rtcm->buff,i,10)*1E-10;        i+=10;
        eph.svh   =getbitu(rtcm->buff,i, 1);              i+= 1;
    }
    else {
        trace(2,"rtcm3 1042 length error: len=%d\n",rtcm->len);
        return -1;
    }
    trace(4,"decode_type1042: prn=%d iode=%d toe=%.0f\n",prn,eph.iode,eph.toes);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," prn=%2d iode=%3d iodc=%3d week=%d toe=%6.0f toc=%6.0f svh=%02X",
                prn,eph.iode,eph.iodc,week,eph.toes,toc,eph.svh);
    }
    if (!(sat=satno(sys,prn))) {
        trace(2,"rtcm3 1042 satellite number error: prn=%d\n",prn);
        return -1;
    }
    eph.sat=sat;
    eph.week=od_rtk_rtcm3_adjbdtweek(week);
    if (rtcm->time.time==0) rtcm->time=utc2gpst(timeget());
    tt=timediff(bdt2gpst(bdt2time(eph.week,eph.toes)),rtcm->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt>=302400.0) eph.week--;
    eph.toe=bdt2gpst(bdt2time(eph.week,eph.toes)); /* bdt -> gpst */
    eph.toc=bdt2gpst(bdt2time(eph.week,toc));      /* bdt -> gpst */
    eph.ttr=rtcm->time;
    eph.A=sqrtA*sqrtA;
    if (!strstr(rtcm->opt,"-EPHALL")) {
        if (timediff(eph.toe,rtcm->nav.eph[sat-1].toe)==0.0&&
            eph.iode==rtcm->nav.eph[sat-1].iode&&
            eph.iodc==rtcm->nav.eph[sat-1].iodc) return 0; /* unchanged */
    }
    rtcm->nav.eph[sat-1]=eph;
    rtcm->ephset=0;
    rtcm->ephsat=sat;
    return 2;
}
/* decode SSR message epoch time ---------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr_epoch(rtcm_t *rtcm, int sys, int subtype)
{
    double tod,tow;
    int i=24+12;
    
    if (subtype==0) { /* RTCM SSR */
        
        if (sys==SYS_GLO) {
            tod=getbitu(rtcm->buff,i,17); i+=17;
            od_rtk_rtcm3_adjday_glot(rtcm,tod);
        }
        else {
            tow=getbitu(rtcm->buff,i,20); i+=20;
            od_rtk_rtcm3_adjweek(rtcm,tow);
        }
    }
    else { /* IGS SSR */
        i+=3+8;
        tow=getbitu(rtcm->buff,i,20); i+=20;
        od_rtk_rtcm3_adjweek(rtcm,tow);
    }
    return i;
}
/* decode SSR 1,4 message header ---------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr1_head(rtcm_t *rtcm, int sys, int subtype, int *sync,
                            int *iod, double *udint, int *refd, int *hsize)
{
    char *msg,tstr[64];
    int i=24+12,nsat,udi,provid=0,solid=0,ns;
    
    if (subtype==0) { /* RTCM SSR */
        ns=(sys==SYS_QZS)?4:6;
        if (i+((sys==SYS_GLO)?53:50+ns)>rtcm->len*8) return -1;
    }
    else { /* IGS SSR */
        ns=6;
        if (i+3+8+50+ns>rtcm->len*8) return -1;
    }
    i=od_rtk_rtcm3_decode_ssr_epoch(rtcm,sys,subtype);
    udi   =getbitu(rtcm->buff,i, 4); i+= 4;
    *sync =getbitu(rtcm->buff,i, 1); i+= 1;
    if (subtype==0) { /* RTCM SSR */
        *refd=getbitu(rtcm->buff,i,1); i+=1; /* satellite ref datum */
    }
    *iod  =getbitu(rtcm->buff,i, 4); i+= 4; /* IOD SSR */
    provid=getbitu(rtcm->buff,i,16); i+=16; /* provider ID */
    solid =getbitu(rtcm->buff,i, 4); i+= 4; /* solution ID */
    if (subtype>0) { /* IGS SSR */
        *refd=getbitu(rtcm->buff,i,1); i+=1; /* global/regional CRS indicator */
    }
    nsat  =getbitu(rtcm->buff,i,ns); i+=ns;
    *udint=od_rtk_rtcm3_ssrudint[udi];
    
    time2str(rtcm->time,tstr,2);
    trace(4,"decode_ssr1_head: time=%s sys=%d subtype=%d nsat=%d sync=%d iod=%d"
         " provid=%d solid=%d\n",tstr,sys,subtype,nsat,*sync,*iod,provid,solid);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," %s nsat=%2d iod=%2d udi=%2d sync=%d",tstr,nsat,*iod,udi,
                *sync);
    }
    *hsize=i;
    return nsat;
}
/* decode SSR 2,3,5,6 message header -----------------------------------------*/
static int od_rtk_rtcm3_decode_ssr2_head(rtcm_t *rtcm, int sys, int subtype, int *sync,
                            int *iod, double *udint, int *hsize)
{
    char *msg,tstr[64];
    int i=24+12,nsat,udi,provid=0,solid=0,ns;
    
    if (subtype==0) { /* RTCM SSR */
        ns=(sys==SYS_QZS)?4:6;
        if (i+((sys==SYS_GLO)?52:49+ns)>rtcm->len*8) return -1;
    }
    else {
        ns=6;
        if (i+3+8+49+ns>rtcm->len*8) return -1;
    }
    i=od_rtk_rtcm3_decode_ssr_epoch(rtcm,sys,subtype);
    udi   =getbitu(rtcm->buff,i, 4); i+= 4;
    *sync =getbitu(rtcm->buff,i, 1); i+= 1;
    *iod  =getbitu(rtcm->buff,i, 4); i+= 4;
    provid=getbitu(rtcm->buff,i,16); i+=16; /* provider ID */
    solid =getbitu(rtcm->buff,i, 4); i+= 4; /* solution ID */
    nsat  =getbitu(rtcm->buff,i,ns); i+=ns;
    *udint=od_rtk_rtcm3_ssrudint[udi];
    
    time2str(rtcm->time,tstr,2);
    trace(4,"decode_ssr2_head: time=%s sys=%d subtype=%d nsat=%d sync=%d iod=%d"
         " provid=%d solid=%d\n",tstr,sys,subtype,nsat,*sync,*iod,provid,solid);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," %s nsat=%2d iod=%2d udi=%2d sync=%d",tstr,nsat,*iod,udi,
                *sync);
    }
    *hsize=i;
    return nsat;
}
/* decode SSR 1: orbit corrections -------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr1(rtcm_t *rtcm, int sys, int subtype)
{
    double udint,deph[3],ddeph[3];
    int i,j,k,type,sync,iod,nsat,prn,sat,iode,iodcrc=0,refd=0,np,ni,nj,offp;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr1_head(rtcm,sys,subtype,&sync,&iod,&udint,&refd,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; ni= 8; nj= 0; offp=  0; break;
        case SYS_GLO: np=5; ni= 8; nj= 0; offp=  0; break;
        case SYS_GAL: np=6; ni=10; nj= 0; offp=  0; break;
        case SYS_QZS: np=4; ni= 8; nj= 0; offp=192; break;
        case SYS_CMP: np=6; ni=10; nj=24; offp=  1; break;
        case SYS_SBS: np=6; ni= 9; nj=24; offp=120; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6; ni=8; nj=0;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+121+np+ni+nj<=rtcm->len*8;j++) {
        prn     =getbitu(rtcm->buff,i,np)+offp; i+=np;
        iode    =getbitu(rtcm->buff,i,ni);      i+=ni;
        iodcrc  =getbitu(rtcm->buff,i,nj);      i+=nj;
        deph [0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        deph [1]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        deph [2]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        ddeph[0]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        ddeph[1]=getbits(rtcm->buff,i,19)*4E-6; i+=19;
        ddeph[2]=getbits(rtcm->buff,i,19)*4E-6; i+=19;
        
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [0]=rtcm->time;
        rtcm->ssr[sat-1].udi[0]=udint;
        rtcm->ssr[sat-1].iod[0]=iod;
        rtcm->ssr[sat-1].iode=iode;     /* SBAS/BDS: toe/t0 modulo */
        rtcm->ssr[sat-1].iodcrc=iodcrc; /* SBAS/BDS: IOD CRC */
        rtcm->ssr[sat-1].refd=refd;
        
        for (k=0;k<3;k++) {
            rtcm->ssr[sat-1].deph [k]=deph [k];
            rtcm->ssr[sat-1].ddeph[k]=ddeph[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}
/* decode SSR 2: clock corrections -------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr2(rtcm_t *rtcm, int sys, int subtype)
{
    double udint,dclk[3];
    int i,j,k,type,sync,iod,nsat,prn,sat,np,offp;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr2_head(rtcm,sys,subtype,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=192; break;
        case SYS_CMP: np=6; offp=  1; break;
        case SYS_SBS: np=6; offp=120; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+70+np<=rtcm->len*8;j++) {
        prn    =getbitu(rtcm->buff,i,np)+offp; i+=np;
        dclk[0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        dclk[1]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        dclk[2]=getbits(rtcm->buff,i,27)*2E-8; i+=27;
        
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [1]=rtcm->time;
        rtcm->ssr[sat-1].udi[1]=udint;
        rtcm->ssr[sat-1].iod[1]=iod;
        
        for (k=0;k<3;k++) {
            rtcm->ssr[sat-1].dclk[k]=dclk[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}
/* decode SSR 3: satellite code biases ---------------------------------------*/
static int od_rtk_rtcm3_decode_ssr3(rtcm_t *rtcm, int sys, int subtype)
{
    const uint8_t *sigs;
    double udint,bias,cbias[MAXCODE];
    int i,j,k,type,mode,sync,iod,nsat,prn,sat,nbias,np,offp;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr2_head(rtcm,sys,subtype,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; sigs=ssr_sig_gps; break;
        case SYS_GLO: np=5; offp=  0; sigs=ssr_sig_glo; break;
        case SYS_GAL: np=6; offp=  0; sigs=ssr_sig_gal; break;
        case SYS_QZS: np=4; offp=192; sigs=ssr_sig_qzs; break;
        case SYS_CMP: np=6; offp=  1; sigs=ssr_sig_cmp; break;
        case SYS_SBS: np=6; offp=120; sigs=ssr_sig_sbs; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+5+np<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i,np)+offp; i+=np;
        nbias=getbitu(rtcm->buff,i, 5);      i+= 5;
        
        for (k=0;k<MAXCODE;k++) cbias[k]=0.0;
        for (k=0;k<nbias&&i+19<=rtcm->len*8;k++) {
            mode=getbitu(rtcm->buff,i, 5);      i+= 5;
            bias=getbits(rtcm->buff,i,14)*0.01; i+=14;
            if (sigs[mode]) {
                cbias[sigs[mode]-1]=(float)bias;
            }
            else {
                trace(2,"rtcm3 %d not supported mode: mode=%d\n",type,mode);
            }
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [4]=rtcm->time;
        rtcm->ssr[sat-1].udi[4]=udint;
        rtcm->ssr[sat-1].iod[4]=iod;
        
        for (k=0;k<MAXCODE;k++) {
            rtcm->ssr[sat-1].cbias[k]=(float)cbias[k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}
/* decode SSR 4: combined orbit and clock corrections ------------------------*/
static int od_rtk_rtcm3_decode_ssr4(rtcm_t *rtcm, int sys, int subtype)
{
    double udint,deph[3],ddeph[3],dclk[3];
    int i,j,k,type,nsat,sync,iod,prn,sat,iode,iodcrc=0,refd=0,np,ni,nj,offp;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr1_head(rtcm,sys,subtype,&sync,&iod,&udint,&refd,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; ni= 8; nj= 0; offp=  0; break;
        case SYS_GLO: np=5; ni= 8; nj= 0; offp=  0; break;
        case SYS_GAL: np=6; ni=10; nj= 0; offp=  0; break;
        case SYS_QZS: np=4; ni= 8; nj= 0; offp=192; break;
        case SYS_CMP: np=6; ni=10; nj=24; offp=  1; break;
        case SYS_SBS: np=6; ni= 9; nj=24; offp=120; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6; ni=8; nj=0;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+191+np+ni+nj<=rtcm->len*8;j++) {
        prn     =getbitu(rtcm->buff,i,np)+offp; i+=np;
        iode    =getbitu(rtcm->buff,i,ni);      i+=ni;
        iodcrc  =getbitu(rtcm->buff,i,nj);      i+=nj;
        deph [0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        deph [1]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        deph [2]=getbits(rtcm->buff,i,20)*4E-4; i+=20;
        ddeph[0]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        ddeph[1]=getbits(rtcm->buff,i,19)*4E-6; i+=19;
        ddeph[2]=getbits(rtcm->buff,i,19)*4E-6; i+=19;
        
        dclk [0]=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        dclk [1]=getbits(rtcm->buff,i,21)*1E-6; i+=21;
        dclk [2]=getbits(rtcm->buff,i,27)*2E-8; i+=27;
        
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [0]=rtcm->ssr[sat-1].t0 [1]=rtcm->time;
        rtcm->ssr[sat-1].udi[0]=rtcm->ssr[sat-1].udi[1]=udint;
        rtcm->ssr[sat-1].iod[0]=rtcm->ssr[sat-1].iod[1]=iod;
        rtcm->ssr[sat-1].iode=iode;
        rtcm->ssr[sat-1].iodcrc=iodcrc;
        rtcm->ssr[sat-1].refd=refd;
        
        for (k=0;k<3;k++) {
            rtcm->ssr[sat-1].deph [k]=deph [k];
            rtcm->ssr[sat-1].ddeph[k]=ddeph[k];
            rtcm->ssr[sat-1].dclk [k]=dclk [k];
        }
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}
/* decode SSR 5: URA ---------------------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr5(rtcm_t *rtcm, int sys, int subtype)
{
    double udint;
    int i,j,type,nsat,sync,iod,prn,sat,ura,np,offp;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr2_head(rtcm,sys,subtype,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=192; break;
        case SYS_CMP: np=6; offp=  1; break;
        case SYS_SBS: np=6; offp=120; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+6+np<=rtcm->len*8;j++) {
        prn=getbitu(rtcm->buff,i,np)+offp; i+=np;
        ura=getbitu(rtcm->buff,i, 6);      i+= 6;
        
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [3]=rtcm->time;
        rtcm->ssr[sat-1].udi[3]=udint;
        rtcm->ssr[sat-1].iod[3]=iod;
        rtcm->ssr[sat-1].ura=ura;
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}
/* decode SSR 6: high rate clock correction ----------------------------------*/
static int od_rtk_rtcm3_decode_ssr6(rtcm_t *rtcm, int sys, int subtype)
{
    double udint,hrclk;
    int i,j,type,nsat,sync,iod,prn,sat,np,offp;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr2_head(rtcm,sys,subtype,&sync,&iod,&udint,&i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=192; break;
        case SYS_CMP: np=6; offp=  1; break;
        case SYS_SBS: np=6; offp=120; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+22+np<=rtcm->len*8;j++) {
        prn  =getbitu(rtcm->buff,i,np)+offp; i+=np;
        hrclk=getbits(rtcm->buff,i,22)*1E-4; i+=22;
        
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [2]=rtcm->time;
        rtcm->ssr[sat-1].udi[2]=udint;
        rtcm->ssr[sat-1].iod[2]=iod;
        rtcm->ssr[sat-1].hrclk=hrclk;
        rtcm->ssr[sat-1].update=1;
    }
    return sync?0:10;
}
/* decode SSR 7 message header -----------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr7_head(rtcm_t *rtcm, int sys, int subtype, int *sync,
                            int *iod, double *udint, int *dispe, int *mw,
                            int *hsize)
{
    char *msg,tstr[64];
    int i=24+12,nsat,udi,provid=0,solid=0,ns;
    
    if (subtype==0) { /* RTCM SSR */
        ns=(sys==SYS_QZS)?4:6;
        if (i+((sys==SYS_GLO)?54:51+ns)>rtcm->len*8) return -1;
    }
    else { /* IGS SSR */
        ns=6;
        if (i+3+8+51+ns>rtcm->len*8) return -1;
    }
    i=od_rtk_rtcm3_decode_ssr_epoch(rtcm,sys,subtype);
    udi   =getbitu(rtcm->buff,i, 4); i+= 4;
    *sync =getbitu(rtcm->buff,i, 1); i+= 1;
    *iod  =getbitu(rtcm->buff,i, 4); i+= 4;
    provid=getbitu(rtcm->buff,i,16); i+=16; /* provider ID */
    solid =getbitu(rtcm->buff,i, 4); i+= 4; /* solution ID */
    *dispe=getbitu(rtcm->buff,i, 1); i+= 1; /* dispersive bias consistency ind */
    *mw   =getbitu(rtcm->buff,i, 1); i+= 1; /* MW consistency indicator */
    nsat  =getbitu(rtcm->buff,i,ns); i+=ns;
    *udint=od_rtk_rtcm3_ssrudint[udi];
    
    time2str(rtcm->time,tstr,2);
    trace(4,"decode_ssr7_head: time=%s sys=%d subtype=%d nsat=%d sync=%d iod=%d"
          " provid=%d solid=%d\n",tstr,sys,subtype,nsat,*sync,*iod,provid,solid);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," %s nsat=%2d iod=%2d udi=%2d sync=%d",tstr,nsat,*iod,udi,
                *sync);
    }
    *hsize=i;
    return nsat;
}
/* decode SSR 7: phase bias --------------------------------------------------*/
static int od_rtk_rtcm3_decode_ssr7(rtcm_t *rtcm, int sys, int subtype)
{
    const uint8_t *sigs;
    double udint,bias,std=0.0,pbias[MAXCODE],stdpb[MAXCODE];
    int i,j,k,type,mode,sync,iod,nsat,prn,sat,nbias,np,mw,offp,sii,swl;
    int dispe,sdc,yaw_ang,yaw_rate;
    
    type=getbitu(rtcm->buff,24,12);
    
    if ((nsat=od_rtk_rtcm3_decode_ssr7_head(rtcm,sys,subtype,&sync,&iod,&udint,&dispe,&mw,
                               &i))<0) {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; sigs=ssr_sig_gps; break;
        case SYS_GLO: np=5; offp=  0; sigs=ssr_sig_glo; break;
        case SYS_GAL: np=6; offp=  0; sigs=ssr_sig_gal; break;
        case SYS_QZS: np=4; offp=192; sigs=ssr_sig_qzs; break;
        case SYS_CMP: np=6; offp=  1; sigs=ssr_sig_cmp; break;
        default: return sync?0:10;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    for (j=0;j<nsat&&i+5+17+np<=rtcm->len*8;j++) {
        prn     =getbitu(rtcm->buff,i,np)+offp; i+=np;
        nbias   =getbitu(rtcm->buff,i, 5);      i+= 5;
        yaw_ang =getbitu(rtcm->buff,i, 9);      i+= 9;
        yaw_rate=getbits(rtcm->buff,i, 8);      i+= 8;
        
        for (k=0;k<MAXCODE;k++) pbias[k]=stdpb[k]=0.0;
        for (k=0;k<nbias&&i+((subtype==0)?49:32)<=rtcm->len*8;k++) {
            mode=getbitu(rtcm->buff,i, 5); i+= 5;
            sii =getbitu(rtcm->buff,i, 1); i+= 1; /* integer-indicator */
            swl =getbitu(rtcm->buff,i, 2); i+= 2; /* WL integer-indicator */
            sdc =getbitu(rtcm->buff,i, 4); i+= 4; /* discontinuity counter */
            bias=getbits(rtcm->buff,i,20); i+=20; /* phase bias (m) */
            if (subtype==0) {
                std=getbitu(rtcm->buff,i,17); i+=17; /* phase bias std-dev (m) */
            }
            if (sigs[mode]) {
                pbias[sigs[mode]-1]=bias*0.0001; /* (m) */
                stdpb[sigs[mode]-1]=std *0.0001; /* (m) */
            }
            else {
                trace(2,"rtcm3 %d not supported mode: mode=%d\n",type,mode);
            }
        }
        if (!(sat=satno(sys,prn))) {
            trace(2,"rtcm3 %d satellite number error: prn=%d\n",type,prn);
            continue;
        }
        rtcm->ssr[sat-1].t0 [5]=rtcm->time;
        rtcm->ssr[sat-1].udi[5]=udint;
        rtcm->ssr[sat-1].iod[5]=iod;
        rtcm->ssr[sat-1].yaw_ang =yaw_ang / 256.0*180.0; /* (deg) */
        rtcm->ssr[sat-1].yaw_rate=yaw_rate/8192.0*180.0; /* (deg/s) */
        
        for (k=0;k<MAXCODE;k++) {
            rtcm->ssr[sat-1].pbias[k]=pbias[k];
            rtcm->ssr[sat-1].stdpb[k]=(float)stdpb[k];
        }
    }
    return 20;
}
/* get signal index ----------------------------------------------------------*/
static void od_rtk_rtcm3_sigindex(int sys, const uint8_t *code, int n, const char *opt,
                     int *idx)
{
    int i,nex,pri,pri_h[8]={0},index[8]={0},ex[32]={0};
    
    /* test code priority */
    for (i=0;i<n;i++) {
        if (!code[i]) continue;
        
        if (idx[i]>=NFREQ) { /* save as extended signal if idx >= NFREQ */
            ex[i]=1;
            continue;
        }
        /* code priority */
        pri=getcodepri(sys,code[i],opt);
        
        /* select highest priority signal */
        if (pri>pri_h[idx[i]]) {
            if (index[idx[i]]) ex[index[idx[i]]-1]=1;
            pri_h[idx[i]]=pri;
            index[idx[i]]=i+1;
        }
        else ex[i]=1;
    }
    /* signal index in obs data */
    for (i=nex=0;i<n;i++) {
        if (ex[i]==0) ;
        else if (nex<NEXOBS) idx[i]=NFREQ+nex++;
        else { /* no space in obs data */
            trace(2,"rtcm msm: no space in obs data sys=%d code=%d\n",sys,code[i]);
            idx[i]=-1;
        }
#if 0 /* for debug */
        trace(2,"sig pos: sys=%d code=%d ex=%d idx=%d\n",sys,code[i],ex[i],idx[i]);
#endif
    }
}
/* save obs data in MSM message ----------------------------------------------*/
static void od_rtk_rtcm3_save_msm_obs(rtcm_t *rtcm, int sys, od_rtk_rtcm3_msm_h_t *h, const double *r,
                         const double *pr, const double *cp, const double *rr,
                         const double *rrf, const double *cnr, const int *lock,
                         const int *ex, const int *half)
{
    const char *sig[32];
    double tt,freq;
    uint8_t code[32];
    char *msm_type="",*q=NULL;
    int i,j,k,type,prn,sat,fcn,index=0,idx[32];
    
    type=getbitu(rtcm->buff,24,12);
    
    switch (sys) {
        case SYS_GPS: msm_type=q=rtcm->msmtype[0]; break;
        case SYS_GLO: msm_type=q=rtcm->msmtype[1]; break;
        case SYS_GAL: msm_type=q=rtcm->msmtype[2]; break;
        case SYS_QZS: msm_type=q=rtcm->msmtype[3]; break;
        case SYS_SBS: msm_type=q=rtcm->msmtype[4]; break;
        case SYS_CMP: msm_type=q=rtcm->msmtype[5]; break;
        case SYS_IRN: msm_type=q=rtcm->msmtype[6]; break;
    }
    /* id to signal */
    for (i=0;i<h->nsig;i++) {
        switch (sys) {
            case SYS_GPS: sig[i]=msm_sig_gps[h->sigs[i]-1]; break;
            case SYS_GLO: sig[i]=msm_sig_glo[h->sigs[i]-1]; break;
            case SYS_GAL: sig[i]=msm_sig_gal[h->sigs[i]-1]; break;
            case SYS_QZS: sig[i]=msm_sig_qzs[h->sigs[i]-1]; break;
            case SYS_SBS: sig[i]=msm_sig_sbs[h->sigs[i]-1]; break;
            case SYS_CMP: sig[i]=msm_sig_cmp[h->sigs[i]-1]; break;
            case SYS_IRN: sig[i]=msm_sig_irn[h->sigs[i]-1]; break;
            default: sig[i]=""; break;
        }
        /* signal to rinex obs type */
        code[i]=obs2code(sig[i]);
        idx[i]=code2idx(sys,code[i]);
        
        if (code[i]!=CODE_NONE) {
            if (q) q+=sprintf(q,"L%s%s",sig[i],i<h->nsig-1?",":"");
        }
        else {
            if (q) q+=sprintf(q,"(%d)%s",h->sigs[i],i<h->nsig-1?",":"");
            
            trace(2,"rtcm3 %d: unknown signal id=%2d\n",type,h->sigs[i]);
        }
    }
    trace(3,"rtcm3 %d: signals=%s\n",type,msm_type);
    
    /* get signal index */
    od_rtk_rtcm3_sigindex(sys,code,h->nsig,rtcm->opt,idx);
    
    for (i=j=0;i<h->nsat;i++) {
        
        prn=h->sats[i];
        if      (sys==SYS_QZS) prn+=MINPRNQZS-1;
        else if (sys==SYS_SBS) prn+=MINPRNSBS-1;
        
        if ((sat=satno(sys,prn))) {
            tt=timediff(rtcm->obs.data[0].time,rtcm->time);
            if (rtcm->obsflag||fabs(tt)>1E-9) {
                rtcm->obs.n=rtcm->obsflag=0;
            }
            index=od_rtk_rtcm3_obsindex(&rtcm->obs,rtcm->time,sat);
        }
        else {
            trace(2,"rtcm3 %d satellite error: prn=%d\n",type,prn);
        }
        fcn=0;
        if (sys==SYS_GLO) {
            fcn=-8; /* no glonass fcn info */
            if (ex&&ex[i]<=13) {
                fcn=ex[i]-7;
                if (!rtcm->nav.glo_fcn[prn-1]) {
                    rtcm->nav.glo_fcn[prn-1]=fcn+8; /* fcn+8 */
                }
            }
            else if (rtcm->nav.geph[prn-1].sat==sat) {
                fcn=rtcm->nav.geph[prn-1].frq;
            }
            else if (rtcm->nav.glo_fcn[prn-1]>0) {
                fcn=rtcm->nav.glo_fcn[prn-1]-8;
            }
        }
        for (k=0;k<h->nsig;k++) {
            if (!h->cellmask[k+i*h->nsig]) continue;
            
            if (sat&&index>=0&&idx[k]>=0) {
                freq=fcn<-7?0.0:code2freq(sys,code[k],fcn);
                
                /* pseudorange (m) */
                if (r[i]!=0.0&&pr[j]>-1E12) {
                    rtcm->obs.data[index].P[idx[k]]=r[i]+pr[j];
                }
                /* carrier-phase (cycle) */
                if (r[i]!=0.0&&cp[j]>-1E12) {
                    rtcm->obs.data[index].L[idx[k]]=(r[i]+cp[j])*freq/CLIGHT;
                }
                /* doppler (hz) */
                if (rr&&rrf&&rrf[j]>-1E12) {
                    rtcm->obs.data[index].D[idx[k]]=
                        (float)(-(rr[i]+rrf[j])*freq/CLIGHT);
                }
                rtcm->obs.data[index].LLI[idx[k]]=
                    od_rtk_rtcm3_lossoflock(rtcm,sat,idx[k],lock[j])+(half[j]?3:0);
                rtcm->obs.data[index].SNR [idx[k]]=(uint16_t)(cnr[j]/SNR_UNIT+0.5);
                rtcm->obs.data[index].code[idx[k]]=code[k];
            }
            j++;
        }
    }
}
/* decode type MSM message header --------------------------------------------*/
static int od_rtk_rtcm3_decode_msm_head(rtcm_t *rtcm, int sys, int *sync, int *iod,
                           od_rtk_rtcm3_msm_h_t *h, int *hsize)
{
    od_rtk_rtcm3_msm_h_t h0={0};
    double tow,tod;
    char *msg,tstr[64];
    int i=24,j,dow,mask,staid,type,ncell=0;
    
    type=getbitu(rtcm->buff,i,12); i+=12;
    
    *h=h0;
    if (i+157<=rtcm->len*8) {
        staid     =getbitu(rtcm->buff,i,12);       i+=12;
        
        if (sys==SYS_GLO) {
            dow   =getbitu(rtcm->buff,i, 3);       i+= 3;
            tod   =getbitu(rtcm->buff,i,27)*0.001; i+=27;
            od_rtk_rtcm3_adjday_glot(rtcm,tod);
        }
        else if (sys==SYS_CMP) {
            tow   =getbitu(rtcm->buff,i,30)*0.001; i+=30;
            tow+=14.0; /* BDT -> GPST */
            od_rtk_rtcm3_adjweek(rtcm,tow);
        }
        else {
            tow   =getbitu(rtcm->buff,i,30)*0.001; i+=30;
            od_rtk_rtcm3_adjweek(rtcm,tow);
        }
        *sync     =getbitu(rtcm->buff,i, 1);       i+= 1;
        *iod      =getbitu(rtcm->buff,i, 3);       i+= 3;
        h->time_s =getbitu(rtcm->buff,i, 7);       i+= 7;
        h->clk_str=getbitu(rtcm->buff,i, 2);       i+= 2;
        h->clk_ext=getbitu(rtcm->buff,i, 2);       i+= 2;
        h->smooth =getbitu(rtcm->buff,i, 1);       i+= 1;
        h->tint_s =getbitu(rtcm->buff,i, 3);       i+= 3;
        for (j=1;j<=64;j++) {
            mask=getbitu(rtcm->buff,i,1); i+=1;
            if (mask) h->sats[h->nsat++]=j;
        }
        for (j=1;j<=32;j++) {
            mask=getbitu(rtcm->buff,i,1); i+=1;
            if (mask) h->sigs[h->nsig++]=j;
        }
    }
    else {
        trace(2,"rtcm3 %d length error: len=%d\n",type,rtcm->len);
        return -1;
    }
    /* test station id */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    if (h->nsat*h->nsig>64) {
        trace(2,"rtcm3 %d number of sats and sigs error: nsat=%d nsig=%d\n",
              type,h->nsat,h->nsig);
        return -1;
    }
    if (i+h->nsat*h->nsig>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: len=%d nsat=%d nsig=%d\n",type,
              rtcm->len,h->nsat,h->nsig);
        return -1;
    }
    for (j=0;j<h->nsat*h->nsig;j++) {
        h->cellmask[j]=getbitu(rtcm->buff,i,1); i+=1;
        if (h->cellmask[j]) ncell++;
    }
    *hsize=i;
    
    time2str(rtcm->time,tstr,2);
    trace(4,"decode_head_msm: time=%s sys=%d staid=%d nsat=%d nsig=%d sync=%d iod=%d ncell=%d\n",
          tstr,sys,staid,h->nsat,h->nsig,*sync,*iod,ncell);
    
    if (rtcm->outtype) {
        msg=rtcm->msgtype+strlen(rtcm->msgtype);
        sprintf(msg," staid=%4d %s nsat=%2d nsig=%2d iod=%2d ncell=%2d sync=%d",
                staid,tstr,h->nsat,h->nsig,*iod,ncell,*sync);
    }
    return ncell;
}
/* decode unsupported MSM message --------------------------------------------*/
static int od_rtk_rtcm3_decode_msm0(rtcm_t *rtcm, int sys)
{
    od_rtk_rtcm3_msm_h_t h={0};
    int i,sync,iod;
    if (od_rtk_rtcm3_decode_msm_head(rtcm,sys,&sync,&iod,&h,&i)<0) return -1;
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode MSM 4: full pseudorange and phaserange plus CNR --------------------*/
static int od_rtk_rtcm3_decode_msm4(rtcm_t *rtcm, int sys)
{
    od_rtk_rtcm3_msm_h_t h={0};
    double r[64],pr[64],cp[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,prv,cpv,lock[64],half[64];
    
    type=getbitu(rtcm->buff,24,12);
    
    /* decode msm header */
    if ((ncell=od_rtk_rtcm3_decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    
    if (i+h.nsat*18+ncell*48>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
              ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) r[j]=0.0;
    for (j=0;j<ncell;j++) pr[j]=cp[j]=-1E16;
    
    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,15); i+=15;
        if (prv!=-16384) pr[j]=prv*P2_24*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,22); i+=22;
        if (cpv!=-2097152) cp[j]=cpv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,4); i+=4;
    }
    for (j=0;j<ncell;j++) { /* half-cycle ambiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,6)*1.0; i+=6;
    }
    /* save obs data in msm message */
    od_rtk_rtcm3_save_msm_obs(rtcm,sys,&h,r,pr,cp,NULL,NULL,cnr,lock,NULL,half);
    
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode MSM 5: full pseudorange, phaserange, phaserangerate and CNR --------*/
static int od_rtk_rtcm3_decode_msm5(rtcm_t *rtcm, int sys)
{
    od_rtk_rtcm3_msm_h_t h={0};
    double r[64],rr[64],pr[64],cp[64],rrf[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,rate,prv,cpv,rrv,lock[64];
    int ex[64],half[64];
    
    type=getbitu(rtcm->buff,24,12);
    
    /* decode msm header */
    if ((ncell=od_rtk_rtcm3_decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    
    if (i+h.nsat*36+ncell*63>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
              ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) {
        r[j]=rr[j]=0.0; ex[j]=15;
    }
    for (j=0;j<ncell;j++) pr[j]=cp[j]=rrf[j]=-1E16;
    
    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* extended info */
        ex[j]=getbitu(rtcm->buff,i, 4); i+= 4;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* phaserangerate */
        rate =getbits(rtcm->buff,i,14); i+=14;
        if (rate!=-8192) rr[j]=rate*1.0;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,15); i+=15;
        if (prv!=-16384) pr[j]=prv*P2_24*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,22); i+=22;
        if (cpv!=-2097152) cp[j]=cpv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,4); i+=4;
    }
    for (j=0;j<ncell;j++) { /* half-cycle ambiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,6)*1.0; i+=6;
    }
    for (j=0;j<ncell;j++) { /* phaserangerate */
        rrv=getbits(rtcm->buff,i,15); i+=15;
        if (rrv!=-16384) rrf[j]=rrv*0.0001;
    }
    /* save obs data in msm message */
    od_rtk_rtcm3_save_msm_obs(rtcm,sys,&h,r,pr,cp,rr,rrf,cnr,lock,ex,half);
    
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode MSM 6: full pseudorange and phaserange plus CNR (high-res) ---------*/
static int od_rtk_rtcm3_decode_msm6(rtcm_t *rtcm, int sys)
{
    od_rtk_rtcm3_msm_h_t h={0};
    double r[64],pr[64],cp[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,prv,cpv,lock[64],half[64];
    
    type=getbitu(rtcm->buff,24,12);
    
    /* decode msm header */
    if ((ncell=od_rtk_rtcm3_decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    
    if (i+h.nsat*18+ncell*65>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
              ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) r[j]=0.0;
    for (j=0;j<ncell;j++) pr[j]=cp[j]=-1E16;
    
    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,20); i+=20;
        if (prv!=-524288) pr[j]=prv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,24); i+=24;
        if (cpv!=-8388608) cp[j]=cpv*P2_31*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,10); i+=10;
    }
    for (j=0;j<ncell;j++) { /* half-cycle ambiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,10)*0.0625; i+=10;
    }
    /* save obs data in msm message */
    od_rtk_rtcm3_save_msm_obs(rtcm,sys,&h,r,pr,cp,NULL,NULL,cnr,lock,NULL,half);
    
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode MSM 7: full pseudorange, phaserange, phaserangerate and CNR (h-res) */
static int od_rtk_rtcm3_decode_msm7(rtcm_t *rtcm, int sys)
{
    od_rtk_rtcm3_msm_h_t h={0};
    double r[64],rr[64],pr[64],cp[64],rrf[64],cnr[64];
    int i,j,type,sync,iod,ncell,rng,rng_m,rate,prv,cpv,rrv,lock[64];
    int ex[64],half[64];
    
    type=getbitu(rtcm->buff,24,12);
    
    /* decode msm header */
    if ((ncell=od_rtk_rtcm3_decode_msm_head(rtcm,sys,&sync,&iod,&h,&i))<0) return -1;
    
    if (i+h.nsat*36+ncell*80>rtcm->len*8) {
        trace(2,"rtcm3 %d length error: nsat=%d ncell=%d len=%d\n",type,h.nsat,
              ncell,rtcm->len);
        return -1;
    }
    for (j=0;j<h.nsat;j++) {
        r[j]=rr[j]=0.0; ex[j]=15;
    }
    for (j=0;j<ncell;j++) pr[j]=cp[j]=rrf[j]=-1E16;
    
    /* decode satellite data */
    for (j=0;j<h.nsat;j++) { /* range */
        rng  =getbitu(rtcm->buff,i, 8); i+= 8;
        if (rng!=255) r[j]=rng*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* extended info */
        ex[j]=getbitu(rtcm->buff,i, 4); i+= 4;
    }
    for (j=0;j<h.nsat;j++) {
        rng_m=getbitu(rtcm->buff,i,10); i+=10;
        if (r[j]!=0.0) r[j]+=rng_m*P2_10*RANGE_MS;
    }
    for (j=0;j<h.nsat;j++) { /* phaserangerate */
        rate =getbits(rtcm->buff,i,14); i+=14;
        if (rate!=-8192) rr[j]=rate*1.0;
    }
    /* decode signal data */
    for (j=0;j<ncell;j++) { /* pseudorange */
        prv=getbits(rtcm->buff,i,20); i+=20;
        if (prv!=-524288) pr[j]=prv*P2_29*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* phaserange */
        cpv=getbits(rtcm->buff,i,24); i+=24;
        if (cpv!=-8388608) cp[j]=cpv*P2_31*RANGE_MS;
    }
    for (j=0;j<ncell;j++) { /* lock time */
        lock[j]=getbitu(rtcm->buff,i,10); i+=10;
    }
    for (j=0;j<ncell;j++) { /* half-cycle amiguity */
        half[j]=getbitu(rtcm->buff,i,1); i+=1;
    }
    for (j=0;j<ncell;j++) { /* cnr */
        cnr[j]=getbitu(rtcm->buff,i,10)*0.0625; i+=10;
    }
    for (j=0;j<ncell;j++) { /* phaserangerate */
        rrv=getbits(rtcm->buff,i,15); i+=15;
        if (rrv!=-16384) rrf[j]=rrv*0.0001;
    }
    /* save obs data in msm message */
    od_rtk_rtcm3_save_msm_obs(rtcm,sys,&h,r,pr,cp,rr,rrf,cnr,lock,ex,half);
    
    rtcm->obsflag=!sync;
    return sync?0:1;
}
/* decode type 1230: GLONASS L1 and L2 code-phase biases ---------------------*/
static int od_rtk_rtcm3_decode_type1230(rtcm_t *rtcm)
{
    int i=24+12,j,staid,align,mask,bias;
    
    if (i+20>=rtcm->len*8) {
        trace(2,"rtcm3 1230: length error len=%d\n",rtcm->len);
        return -1;
    }
    staid=getbitu(rtcm->buff,i,12); i+=12;
    align=getbitu(rtcm->buff,i, 1); i+= 1+3;
    mask =getbitu(rtcm->buff,i, 4); i+= 4;
    
    if (rtcm->outtype) {
        sprintf(rtcm->msgtype+strlen(rtcm->msgtype),
                " staid=%4d align=%d mask=0x%X",staid,align,mask);
    }
    /* test station ID */
    if (!od_rtk_rtcm3_test_staid(rtcm,staid)) return -1;
    
    rtcm->sta.glo_cp_align=align;
    for (j=0;j<4;j++) {
        rtcm->sta.glo_cp_bias[j]=0.0;
    }
    for (j=0;j<4&&i+16<=rtcm->len*8;j++) {
        if (!(mask&(1<<(3-j)))) continue;
        bias=getbits(rtcm->buff,i,16); i+=16;
        if (bias!=-32768) {
            rtcm->sta.glo_cp_bias[j]=bias*0.02;
        }
    }
    return 5;
}
/* decode type 4073: proprietary message Mitsubishi Electric -----------------*/
static int od_rtk_rtcm3_decode_type4073(rtcm_t *rtcm)
{
    int i=24+12,subtype;
    
    subtype=getbitu(rtcm->buff,i,4); i+=4;
    
    if (rtcm->outtype) {
        sprintf(rtcm->msgtype+strlen(rtcm->msgtype)," subtype=%d",subtype);
    }
    trace(2,"rtcm3 4073: unsupported message subtype=%d\n",subtype);
    return 0;
}
/* decode type 4076: proprietary message IGS ---------------------------------*/
static int od_rtk_rtcm3_decode_type4076(rtcm_t *rtcm)
{
    int i=24+12,ver,subtype;
    
    if (i+3+8>=rtcm->len*8) {
        trace(2,"rtcm3 4076: length error len=%d\n",rtcm->len);
        return -1;
    }
    ver    =getbitu(rtcm->buff,i,3); i+=3;
    subtype=getbitu(rtcm->buff,i,8); i+=8;
    
    if (rtcm->outtype) {
        sprintf(rtcm->msgtype+strlen(rtcm->msgtype)," ver=%d subtype=%3d",ver,
                subtype);
    }
    switch (subtype) {
        case  21: return od_rtk_rtcm3_decode_ssr1(rtcm,SYS_GPS,subtype);
        case  22: return od_rtk_rtcm3_decode_ssr2(rtcm,SYS_GPS,subtype);
        case  23: return od_rtk_rtcm3_decode_ssr4(rtcm,SYS_GPS,subtype);
        case  24: return od_rtk_rtcm3_decode_ssr6(rtcm,SYS_GPS,subtype);
        case  25: return od_rtk_rtcm3_decode_ssr3(rtcm,SYS_GPS,subtype);
        case  26: return od_rtk_rtcm3_decode_ssr7(rtcm,SYS_GPS,subtype);
        case  27: return od_rtk_rtcm3_decode_ssr5(rtcm,SYS_GPS,subtype);
        case  41: return od_rtk_rtcm3_decode_ssr1(rtcm,SYS_GLO,subtype);
        case  42: return od_rtk_rtcm3_decode_ssr2(rtcm,SYS_GLO,subtype);
        case  43: return od_rtk_rtcm3_decode_ssr4(rtcm,SYS_GLO,subtype);
        case  44: return od_rtk_rtcm3_decode_ssr6(rtcm,SYS_GLO,subtype);
        case  45: return od_rtk_rtcm3_decode_ssr3(rtcm,SYS_GLO,subtype);
        case  46: return od_rtk_rtcm3_decode_ssr7(rtcm,SYS_GLO,subtype);
        case  47: return od_rtk_rtcm3_decode_ssr5(rtcm,SYS_GLO,subtype);
        case  61: return od_rtk_rtcm3_decode_ssr1(rtcm,SYS_GAL,subtype);
        case  62: return od_rtk_rtcm3_decode_ssr2(rtcm,SYS_GAL,subtype);
        case  63: return od_rtk_rtcm3_decode_ssr4(rtcm,SYS_GAL,subtype);
        case  64: return od_rtk_rtcm3_decode_ssr6(rtcm,SYS_GAL,subtype);
        case  65: return od_rtk_rtcm3_decode_ssr3(rtcm,SYS_GAL,subtype);
        case  66: return od_rtk_rtcm3_decode_ssr7(rtcm,SYS_GAL,subtype);
        case  67: return od_rtk_rtcm3_decode_ssr5(rtcm,SYS_GAL,subtype);
        case  81: return od_rtk_rtcm3_decode_ssr1(rtcm,SYS_QZS,subtype);
        case  82: return od_rtk_rtcm3_decode_ssr2(rtcm,SYS_QZS,subtype);
        case  83: return od_rtk_rtcm3_decode_ssr4(rtcm,SYS_QZS,subtype);
        case  84: return od_rtk_rtcm3_decode_ssr6(rtcm,SYS_QZS,subtype);
        case  85: return od_rtk_rtcm3_decode_ssr3(rtcm,SYS_QZS,subtype);
        case  86: return od_rtk_rtcm3_decode_ssr7(rtcm,SYS_QZS,subtype);
        case  87: return od_rtk_rtcm3_decode_ssr5(rtcm,SYS_QZS,subtype);
        case 101: return od_rtk_rtcm3_decode_ssr1(rtcm,SYS_CMP,subtype);
        case 102: return od_rtk_rtcm3_decode_ssr2(rtcm,SYS_CMP,subtype);
        case 103: return od_rtk_rtcm3_decode_ssr4(rtcm,SYS_CMP,subtype);
        case 104: return od_rtk_rtcm3_decode_ssr6(rtcm,SYS_CMP,subtype);
        case 105: return od_rtk_rtcm3_decode_ssr3(rtcm,SYS_CMP,subtype);
        case 106: return od_rtk_rtcm3_decode_ssr7(rtcm,SYS_CMP,subtype);
        case 107: return od_rtk_rtcm3_decode_ssr5(rtcm,SYS_CMP,subtype);
        case 121: return od_rtk_rtcm3_decode_ssr1(rtcm,SYS_SBS,subtype);
        case 122: return od_rtk_rtcm3_decode_ssr2(rtcm,SYS_SBS,subtype);
        case 123: return od_rtk_rtcm3_decode_ssr4(rtcm,SYS_SBS,subtype);
        case 124: return od_rtk_rtcm3_decode_ssr6(rtcm,SYS_SBS,subtype);
        case 125: return od_rtk_rtcm3_decode_ssr3(rtcm,SYS_SBS,subtype);
        case 126: return od_rtk_rtcm3_decode_ssr7(rtcm,SYS_SBS,subtype);
        case 127: return od_rtk_rtcm3_decode_ssr5(rtcm,SYS_SBS,subtype);
    }
    trace(2,"rtcm3 4076: unsupported message subtype=%d\n",subtype);
    return 0;
}
/* decode RTCM ver.3 message -------------------------------------------------*/
extern int decode_rtcm3(rtcm_t *rtcm)
{
    double tow;
    int ret=0,type=getbitu(rtcm->buff,24,12),week;
    
    trace(3,"decode_rtcm3: len=%3d type=%d\n",rtcm->len,type);
    
    if (rtcm->outtype) {
        sprintf(rtcm->msgtype,"RTCM %4d (%4d):",type,rtcm->len);
    }
    /* real-time input option */
    if (strstr(rtcm->opt,"-RT_INP")) {
        tow=time2gpst(utc2gpst(timeget()),&week);
        rtcm->time=gpst2time(week,floor(tow));
    }
    switch (type) {
        case 1001: ret=od_rtk_rtcm3_decode_type1001(rtcm); break; /* not supported */
        case 1002: ret=od_rtk_rtcm3_decode_type1002(rtcm); break;
        case 1003: ret=od_rtk_rtcm3_decode_type1003(rtcm); break; /* not supported */
        case 1004: ret=od_rtk_rtcm3_decode_type1004(rtcm); break;
        case 1005: ret=od_rtk_rtcm3_decode_type1005(rtcm); break;
        case 1006: ret=od_rtk_rtcm3_decode_type1006(rtcm); break;
        case 1007: ret=od_rtk_rtcm3_decode_type1007(rtcm); break;
        case 1008: ret=od_rtk_rtcm3_decode_type1008(rtcm); break;
        case 1009: ret=od_rtk_rtcm3_decode_type1009(rtcm); break; /* not supported */
        case 1010: ret=od_rtk_rtcm3_decode_type1010(rtcm); break;
        case 1011: ret=od_rtk_rtcm3_decode_type1011(rtcm); break; /* not supported */
        case 1012: ret=od_rtk_rtcm3_decode_type1012(rtcm); break;
        case 1013: ret=od_rtk_rtcm3_decode_type1013(rtcm); break; /* not supported */
        case 1019: ret=od_rtk_rtcm3_decode_type1019(rtcm); break;
        case 1020: ret=od_rtk_rtcm3_decode_type1020(rtcm); break;
        case 1021: ret=od_rtk_rtcm3_decode_type1021(rtcm); break; /* not supported */
        case 1022: ret=od_rtk_rtcm3_decode_type1022(rtcm); break; /* not supported */
        case 1023: ret=od_rtk_rtcm3_decode_type1023(rtcm); break; /* not supported */
        case 1024: ret=od_rtk_rtcm3_decode_type1024(rtcm); break; /* not supported */
        case 1025: ret=od_rtk_rtcm3_decode_type1025(rtcm); break; /* not supported */
        case 1026: ret=od_rtk_rtcm3_decode_type1026(rtcm); break; /* not supported */
        case 1027: ret=od_rtk_rtcm3_decode_type1027(rtcm); break; /* not supported */
        case 1029: ret=od_rtk_rtcm3_decode_type1029(rtcm); break;
        case 1030: ret=od_rtk_rtcm3_decode_type1030(rtcm); break; /* not supported */
        case 1031: ret=od_rtk_rtcm3_decode_type1031(rtcm); break; /* not supported */
        case 1032: ret=od_rtk_rtcm3_decode_type1032(rtcm); break; /* not supported */
        case 1033: ret=od_rtk_rtcm3_decode_type1033(rtcm); break;
        case 1034: ret=od_rtk_rtcm3_decode_type1034(rtcm); break; /* not supported */
        case 1035: ret=od_rtk_rtcm3_decode_type1035(rtcm); break; /* not supported */
        case 1037: ret=od_rtk_rtcm3_decode_type1037(rtcm); break; /* not supported */
        case 1038: ret=od_rtk_rtcm3_decode_type1038(rtcm); break; /* not supported */
        case 1039: ret=od_rtk_rtcm3_decode_type1039(rtcm); break; /* not supported */
        case 1041: ret=od_rtk_rtcm3_decode_type1041(rtcm); break;
        case 1044: ret=od_rtk_rtcm3_decode_type1044(rtcm); break;
        case 1045: ret=od_rtk_rtcm3_decode_type1045(rtcm); break;
        case 1046: ret=od_rtk_rtcm3_decode_type1046(rtcm); break;
        case   63: ret=od_rtk_rtcm3_decode_type1042(rtcm); break; /* RTCM draft */
        case 1042: ret=od_rtk_rtcm3_decode_type1042(rtcm); break;
        case 1057: ret=od_rtk_rtcm3_decode_ssr1(rtcm,SYS_GPS,0); break;
        case 1058: ret=od_rtk_rtcm3_decode_ssr2(rtcm,SYS_GPS,0); break;
        case 1059: ret=od_rtk_rtcm3_decode_ssr3(rtcm,SYS_GPS,0); break;
        case 1060: ret=od_rtk_rtcm3_decode_ssr4(rtcm,SYS_GPS,0); break;
        case 1061: ret=od_rtk_rtcm3_decode_ssr5(rtcm,SYS_GPS,0); break;
        case 1062: ret=od_rtk_rtcm3_decode_ssr6(rtcm,SYS_GPS,0); break;
        case 1063: ret=od_rtk_rtcm3_decode_ssr1(rtcm,SYS_GLO,0); break;
        case 1064: ret=od_rtk_rtcm3_decode_ssr2(rtcm,SYS_GLO,0); break;
        case 1065: ret=od_rtk_rtcm3_decode_ssr3(rtcm,SYS_GLO,0); break;
        case 1066: ret=od_rtk_rtcm3_decode_ssr4(rtcm,SYS_GLO,0); break;
        case 1067: ret=od_rtk_rtcm3_decode_ssr5(rtcm,SYS_GLO,0); break;
        case 1068: ret=od_rtk_rtcm3_decode_ssr6(rtcm,SYS_GLO,0); break;
        case 1071: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GPS); break; /* not supported */
        case 1072: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GPS); break; /* not supported */
        case 1073: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GPS); break; /* not supported */
        case 1074: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_GPS); break;
        case 1075: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_GPS); break;
        case 1076: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_GPS); break;
        case 1077: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_GPS); break;
        case 1081: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GLO); break; /* not supported */
        case 1082: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GLO); break; /* not supported */
        case 1083: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GLO); break; /* not supported */
        case 1084: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_GLO); break;
        case 1085: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_GLO); break;
        case 1086: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_GLO); break;
        case 1087: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_GLO); break;
        case 1091: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GAL); break; /* not supported */
        case 1092: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GAL); break; /* not supported */
        case 1093: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_GAL); break; /* not supported */
        case 1094: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_GAL); break;
        case 1095: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_GAL); break;
        case 1096: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_GAL); break;
        case 1097: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_GAL); break;
        case 1101: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_SBS); break; /* not supported */
        case 1102: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_SBS); break; /* not supported */
        case 1103: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_SBS); break; /* not supported */
        case 1104: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_SBS); break;
        case 1105: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_SBS); break;
        case 1106: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_SBS); break;
        case 1107: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_SBS); break;
        case 1111: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_QZS); break; /* not supported */
        case 1112: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_QZS); break; /* not supported */
        case 1113: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_QZS); break; /* not supported */
        case 1114: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_QZS); break;
        case 1115: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_QZS); break;
        case 1116: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_QZS); break;
        case 1117: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_QZS); break;
        case 1121: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_CMP); break; /* not supported */
        case 1122: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_CMP); break; /* not supported */
        case 1123: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_CMP); break; /* not supported */
        case 1124: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_CMP); break;
        case 1125: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_CMP); break;
        case 1126: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_CMP); break;
        case 1127: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_CMP); break;
        case 1131: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_IRN); break; /* not supported */
        case 1132: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_IRN); break; /* not supported */
        case 1133: ret=od_rtk_rtcm3_decode_msm0(rtcm,SYS_IRN); break; /* not supported */
        case 1134: ret=od_rtk_rtcm3_decode_msm4(rtcm,SYS_IRN); break;
        case 1135: ret=od_rtk_rtcm3_decode_msm5(rtcm,SYS_IRN); break;
        case 1136: ret=od_rtk_rtcm3_decode_msm6(rtcm,SYS_IRN); break;
        case 1137: ret=od_rtk_rtcm3_decode_msm7(rtcm,SYS_IRN); break;
        case 1230: ret=od_rtk_rtcm3_decode_type1230(rtcm);     break;
        case 1240: ret=od_rtk_rtcm3_decode_ssr1(rtcm,SYS_GAL,0); break; /* draft */
        case 1241: ret=od_rtk_rtcm3_decode_ssr2(rtcm,SYS_GAL,0); break; /* draft */
        case 1242: ret=od_rtk_rtcm3_decode_ssr3(rtcm,SYS_GAL,0); break; /* draft */
        case 1243: ret=od_rtk_rtcm3_decode_ssr4(rtcm,SYS_GAL,0); break; /* draft */
        case 1244: ret=od_rtk_rtcm3_decode_ssr5(rtcm,SYS_GAL,0); break; /* draft */
        case 1245: ret=od_rtk_rtcm3_decode_ssr6(rtcm,SYS_GAL,0); break; /* draft */
        case 1246: ret=od_rtk_rtcm3_decode_ssr1(rtcm,SYS_QZS,0); break; /* draft */
        case 1247: ret=od_rtk_rtcm3_decode_ssr2(rtcm,SYS_QZS,0); break; /* draft */
        case 1248: ret=od_rtk_rtcm3_decode_ssr3(rtcm,SYS_QZS,0); break; /* draft */
        case 1249: ret=od_rtk_rtcm3_decode_ssr4(rtcm,SYS_QZS,0); break; /* draft */
        case 1250: ret=od_rtk_rtcm3_decode_ssr5(rtcm,SYS_QZS,0); break; /* draft */
        case 1251: ret=od_rtk_rtcm3_decode_ssr6(rtcm,SYS_QZS,0); break; /* draft */
        case 1252: ret=od_rtk_rtcm3_decode_ssr1(rtcm,SYS_SBS,0); break; /* draft */
        case 1253: ret=od_rtk_rtcm3_decode_ssr2(rtcm,SYS_SBS,0); break; /* draft */
        case 1254: ret=od_rtk_rtcm3_decode_ssr3(rtcm,SYS_SBS,0); break; /* draft */
        case 1255: ret=od_rtk_rtcm3_decode_ssr4(rtcm,SYS_SBS,0); break; /* draft */
        case 1256: ret=od_rtk_rtcm3_decode_ssr5(rtcm,SYS_SBS,0); break; /* draft */
        case 1257: ret=od_rtk_rtcm3_decode_ssr6(rtcm,SYS_SBS,0); break; /* draft */
        case 1258: ret=od_rtk_rtcm3_decode_ssr1(rtcm,SYS_CMP,0); break; /* draft */
        case 1259: ret=od_rtk_rtcm3_decode_ssr2(rtcm,SYS_CMP,0); break; /* draft */
        case 1260: ret=od_rtk_rtcm3_decode_ssr3(rtcm,SYS_CMP,0); break; /* draft */
        case 1261: ret=od_rtk_rtcm3_decode_ssr4(rtcm,SYS_CMP,0); break; /* draft */
        case 1262: ret=od_rtk_rtcm3_decode_ssr5(rtcm,SYS_CMP,0); break; /* draft */
        case 1263: ret=od_rtk_rtcm3_decode_ssr6(rtcm,SYS_CMP,0); break; /* draft */
        case   11: ret=od_rtk_rtcm3_decode_ssr7(rtcm,SYS_GPS,0); break; /* tentative */
        case   12: ret=od_rtk_rtcm3_decode_ssr7(rtcm,SYS_GAL,0); break; /* tentative */
        case   13: ret=od_rtk_rtcm3_decode_ssr7(rtcm,SYS_QZS,0); break; /* tentative */
        case   14: ret=od_rtk_rtcm3_decode_ssr7(rtcm,SYS_CMP,0); break; /* tentative */
        case 4073: ret=od_rtk_rtcm3_decode_type4073(rtcm); break;
        case 4076: ret=od_rtk_rtcm3_decode_type4076(rtcm); break;
    }
    if (ret>=0) {
        if      (1001<=type&&type<=1299) rtcm->nmsg3[type-1000]++; /*   1-299 */
        else if (4070<=type&&type<=4099) rtcm->nmsg3[type-3770]++; /* 300-329 */
        else rtcm->nmsg3[0]++; /* other */
    }
    return ret;
}


#pragma pop_macro("RANGE_MS")
#pragma pop_macro("PRUNIT_GPS")
#pragma pop_macro("PRUNIT_GLO")
#pragma pop_macro("P2_66")
#pragma pop_macro("P2_59")
#pragma pop_macro("P2_46")
#pragma pop_macro("P2_41")
#pragma pop_macro("P2_34")
#pragma pop_macro("P2_28")
#pragma pop_macro("P2_10")


/* ===== Embedded rtcm3e.c ===== */
#pragma push_macro("MIN")
#undef MIN
#pragma push_macro("P2_10")
#undef P2_10
#pragma push_macro("P2_28")
#undef P2_28
#pragma push_macro("P2_34")
#undef P2_34
#pragma push_macro("P2_41")
#undef P2_41
#pragma push_macro("P2_46")
#undef P2_46
#pragma push_macro("P2_59")
#undef P2_59
#pragma push_macro("P2_66")
#undef P2_66
#pragma push_macro("PRUNIT_GLO")
#undef PRUNIT_GLO
#pragma push_macro("PRUNIT_GPS")
#undef PRUNIT_GPS
#pragma push_macro("RANGE_MS")
#undef RANGE_MS
#pragma push_macro("ROUND")
#undef ROUND
#pragma push_macro("ROUND_U")
#undef ROUND_U

/*------------------------------------------------------------------------------
* rtcm3e.c : rtcm ver.3 message encoder functions
*
*          Copyright (C) 2012-2020 by T.TAKASU, All rights reserved.
*
* references :
*     see rtcm.c
*
* version : $Revision:$ $Date:$
* history : 2012/12/05 1.0  new
*           2012/12/16 1.1  fix bug on ssr high rate clock correction
*           2012/12/24 1.2  fix bug on msm carrier-phase offset correction
*                           fix bug on SBAS sat id in 1001-1004
*                           fix bug on carrier-phase in 1001-1004,1009-1012
*           2012/12/28 1.3  fix bug on compass carrier wave length
*           2013/01/18 1.4  fix bug on ssr message generation
*           2013/05/11 1.5  change type of arg value of setbig()
*           2013/05/19 1.5  gpst -> bdt of time-tag in beidou msm message
*           2013/04/27 1.7  comply with rtcm 3.2 with amendment 1/2 (ref[15])
*                           delete MT 1046 according to ref [15]
*           2014/05/15 1.8  set NT field in MT 1020 glonass ephemeris
*           2014/12/06 1.9  support SBAS/BeiDou SSR messages (ref [16])
*                           fix bug on invalid staid in qzss ssr messages
*           2015/03/22 1.9  add handling of iodcrc for beidou/sbas ssr messages
*           2015/08/03 1.10 fix bug on wrong udint and iod in ssr 7.
*                           support rtcm ssr fcb message mt 2065-2069.
*           2015/09/07 1.11 add message count of MT 2000-2099
*           2015/10/21 1.12 add MT1046 support for IGS MGEX
*           2015/12/04 1.13 add MT63 beidou ephemeris (rtcm draft)
*                           fix bug on msm message generation of beidou
*                           fix bug on ssr 3 message generation (#321)
*           2016/06/12 1.14 fix bug on segmentation fault by generating msm1
*           2016/09/20 1.15 fix bug on MT1045 Galileo week rollover
*           2017/04/11 1.16 fix bug on gst-week in MT1045/1046
*           2018/10/10 1.17 merge changes for 2.4.2 p13
*                           change mt for ssr 7 phase biases
*           2019/05/10 1.21 save galileo E5b data to obs index 2
*           2020/11/30 1.22 support MT1230 GLONASS code-phase biases
*                           support MT1131-1137,1041 (NavIC MSM and ephemeris)
*                           support MT4076 IGS SSR
*                           fixed invalid delta clock C2 value for SSR 2 and 4
*                           delete SSR signal and tracking mode ID table
*                           use API code2idx() to get freq-index
*                           use API code2freq() to get carrier frequency
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/

/* constants and macros ------------------------------------------------------*/

#define PRUNIT_GPS  299792.458          /* rtcm 3 unit of gps pseudorange (m) */
#define PRUNIT_GLO  599584.916          /* rtcm 3 unit of glo pseudorange (m) */
#define RANGE_MS    (CLIGHT*0.001)      /* range in 1 ms */
#define P2_10       0.0009765625          /* 2^-10 */
#define P2_28       3.725290298461914E-09 /* 2^-28 */
#define P2_34       5.820766091346740E-11 /* 2^-34 */
#define P2_41       4.547473508864641E-13 /* 2^-41 */
#define P2_46       1.421085471520200E-14 /* 2^-46 */
#define P2_59       1.734723475976810E-18 /* 2^-59 */
#define P2_66       1.355252715606880E-20 /* 2^-66 */

#define ROUND(x)    ((int)floor((x)+0.5))
#define ROUND_U(x)  ((uint32_t)floor((x)+0.5))
#define MIN(x,y)    ((x)<(y)?(x):(y))

/* MSM signal ID table -------------------------------------------------------*/
extern const char *msm_sig_gps[32];
extern const char *msm_sig_glo[32];
extern const char *msm_sig_gal[32];
extern const char *msm_sig_qzs[32];
extern const char *msm_sig_sbs[32];
extern const char *msm_sig_cmp[32];
extern const char *msm_sig_irn[32];

/* SSR signal and tracking mode IDs ------------------------------------------*/
extern const uint8_t ssr_sig_gps[32];
extern const uint8_t ssr_sig_glo[32];
extern const uint8_t ssr_sig_gal[32];
extern const uint8_t ssr_sig_qzs[32];
extern const uint8_t ssr_sig_cmp[32];
extern const uint8_t ssr_sig_sbs[32];

/* SSR update intervals ------------------------------------------------------*/
static const double od_rtk_rtcm3e_ssrudint[16]={
    1,2,5,10,15,30,60,120,240,300,600,900,1800,3600,7200,10800
};
/* set sign-magnitude bits ---------------------------------------------------*/
static void od_rtk_rtcm3e_setbitg(uint8_t *buff, int pos, int len, int32_t value)
{
    setbitu(buff,pos,1,value<0?1:0);
    setbitu(buff,pos+1,len-1,value<0?-value:value);
}
/* set signed 38 bit field ---------------------------------------------------*/
static void od_rtk_rtcm3e_set38bits(uint8_t *buff, int pos, double value)
{
    int word_h=(int)floor(value/64.0);
    uint32_t word_l=(uint32_t)(value-word_h*64.0);
    setbits(buff,pos  ,32,word_h);
    setbitu(buff,pos+32,6,word_l);
}
/* lock time -----------------------------------------------------------------*/
static int od_rtk_rtcm3e_locktime(gtime_t time, gtime_t *lltime, uint8_t LLI)
{
    if (!lltime->time||(LLI&1)) *lltime=time;
    return (int)timediff(time,*lltime);
}
/* lock time in double -------------------------------------------------------*/
static double od_rtk_rtcm3e_locktime_d(gtime_t time, gtime_t *lltime, uint8_t LLI)
{
    if (!lltime->time||(LLI&1)) *lltime=time;
    return timediff(time,*lltime);
}
/* GLONASS frequency channel number in RTCM (FCN+7,-1:error) -----------------*/
static int od_rtk_rtcm3e_fcn_glo(int sat, rtcm_t *rtcm)
{
    int prn;
    
    if (satsys(sat,&prn)!=SYS_GLO) {
        return -1;
    }
    if (rtcm->nav.geph[prn-1].sat==sat) {
        return rtcm->nav.geph[prn-1].frq+7;
    }
    if (rtcm->nav.glo_fcn[prn-1]>0) { /* fcn+8 (0: no data) */
        return rtcm->nav.glo_fcn[prn-1]-8+7;
    }
    return -1;
}
/* lock time indicator (ref [17] table 3.4-2) --------------------------------*/
static int od_rtk_rtcm3e_to_lock(int lock)
{
    if (lock<0  ) return 0;
    if (lock<24 ) return lock;
    if (lock<72 ) return (lock+24  )/2;
    if (lock<168) return (lock+120 )/4;
    if (lock<360) return (lock+408 )/8;
    if (lock<744) return (lock+1176)/16;
    if (lock<937) return (lock+3096)/32;
    return 127;
}
/* MSM lock time indicator (ref [17] table 3.5-74) ---------------------------*/
static int od_rtk_rtcm3e_to_msm_lock(double lock)
{
    if (lock<0.032  ) return 0;
    if (lock<0.064  ) return 1;
    if (lock<0.128  ) return 2;
    if (lock<0.256  ) return 3;
    if (lock<0.512  ) return 4;
    if (lock<1.024  ) return 5;
    if (lock<2.048  ) return 6;
    if (lock<4.096  ) return 7;
    if (lock<8.192  ) return 8;
    if (lock<16.384 ) return 9;
    if (lock<32.768 ) return 10;
    if (lock<65.536 ) return 11;
    if (lock<131.072) return 12;
    if (lock<262.144) return 13;
    if (lock<524.288) return 14;
    return 15;
}
/* MSM lock time indicator with extended-resolution (ref [17] table 3.5-76) --*/
static int od_rtk_rtcm3e_to_msm_lock_ex(double lock)
{
    int lock_ms = (int)(lock * 1000.0);
    
    if (lock<0.0      ) return 0;
    if (lock<0.064    ) return lock_ms;
    if (lock<0.128    ) return (lock_ms+64       )/2;
    if (lock<0.256    ) return (lock_ms+256      )/4;
    if (lock<0.512    ) return (lock_ms+768      )/8;
    if (lock<1.024    ) return (lock_ms+2048     )/16;
    if (lock<2.048    ) return (lock_ms+5120     )/32;
    if (lock<4.096    ) return (lock_ms+12288    )/64;
    if (lock<8.192    ) return (lock_ms+28672    )/128;
    if (lock<16.384   ) return (lock_ms+65536    )/256;
    if (lock<32.768   ) return (lock_ms+147456   )/512;
    if (lock<65.536   ) return (lock_ms+327680   )/1024;
    if (lock<131.072  ) return (lock_ms+720896   )/2048;
    if (lock<262.144  ) return (lock_ms+1572864  )/4096;
    if (lock<524.288  ) return (lock_ms+3407872  )/8192;
    if (lock<1048.576 ) return (lock_ms+7340032  )/16384;
    if (lock<2097.152 ) return (lock_ms+15728640 )/32768;
    if (lock<4194.304 ) return (lock_ms+33554432 )/65536;
    if (lock<8388.608 ) return (lock_ms+71303168 )/131072;
    if (lock<16777.216) return (lock_ms+150994944)/262144;
    if (lock<33554.432) return (lock_ms+318767104)/524288;
    if (lock<67108.864) return (lock_ms+671088640)/1048576;
    return 704;
}
/* L1 code indicator GPS -----------------------------------------------------*/
static int od_rtk_rtcm3e_to_code1_gps(uint8_t code)
{
    switch (code) {
        case CODE_L1C: return 0; /* L1 C/A */
        case CODE_L1P:
        case CODE_L1W:
        case CODE_L1Y:
        case CODE_L1N: return 1; /* L1 P(Y) direct */
    }
    return 0;
}
/* L2 code indicator GPS -----------------------------------------------------*/
static int od_rtk_rtcm3e_to_code2_gps(uint8_t code)
{
    switch (code) {
        case CODE_L2C:
        case CODE_L2S:
        case CODE_L2L:
        case CODE_L2X: return 0; /* L2 C/A or L2C */
        case CODE_L2P:
        case CODE_L2Y: return 1; /* L2 P(Y) direct */
        case CODE_L2D: return 2; /* L2 P(Y) cross-correlated */
        case CODE_L2W:
        case CODE_L2N: return 3; /* L2 correlated P/Y */
    }
    return 0;
}
/* L1 code indicator GLONASS -------------------------------------------------*/
static int od_rtk_rtcm3e_to_code1_glo(uint8_t code)
{
    switch (code) {
        case CODE_L1C: return 0; /* L1 C/A */
        case CODE_L1P: return 1; /* L1 P */
    }
    return 0;
}
/* L2 code indicator GLONASS -------------------------------------------------*/
static int od_rtk_rtcm3e_to_code2_glo(uint8_t code)
{
    switch (code) {
        case CODE_L2C: return 0; /* L2 C/A */
        case CODE_L2P: return 1; /* L2 P */
    }
    return 0;
}
/* carrier-phase - pseudorange in cycle --------------------------------------*/
static double od_rtk_rtcm3e_cp_pr(double cp, double pr_cyc)
{
    return fmod(cp-pr_cyc+750.0,1500.0)-750.0;
}
/* generate obs field data GPS -----------------------------------------------*/
static void od_rtk_rtcm3e_gen_obs_gps(rtcm_t *rtcm, const obsd_t *data, int *code1, int *pr1,
                        int *ppr1, int *lock1, int *amb, int *cnr1, int *code2,
                        int *pr21, int *ppr2, int *lock2, int *cnr2)
{
    double lam1,lam2,pr1c=0.0,ppr;
    int lt1,lt2;
    
    lam1=CLIGHT/FREQ1;
    lam2=CLIGHT/FREQ2;
    *pr1=*amb=0;
    if (ppr1) *ppr1=0xFFF80000; /* invalid values */
    if (pr21) *pr21=0xFFFFE000;
    if (ppr2) *ppr2=0xFFF80000;
    
    /* L1 peudorange */
    if (data->P[0]!=0.0&&data->code[0]) {
        *amb=(int)floor(data->P[0]/PRUNIT_GPS);
        *pr1=ROUND((data->P[0]-*amb*PRUNIT_GPS)/0.02);
        pr1c=*pr1*0.02+*amb*PRUNIT_GPS;
    }
    /* L1 phaserange - L1 pseudorange */
    if (data->P[0]!=0.0&&data->L[0]!=0.0&&data->code[0]) {
        ppr=od_rtk_rtcm3e_cp_pr(data->L[0],pr1c/lam1);
        if (ppr1) *ppr1=ROUND(ppr*lam1/0.0005);
    }
    /* L2 -L1 pseudorange */
    if (data->P[0]!=0.0&&data->P[1]!=0.0&&data->code[0]&&data->code[1]&&
        fabs(data->P[1]-pr1c)<=163.82) {
        if (pr21) *pr21=ROUND((data->P[1]-pr1c)/0.02);
    }
    /* L2 phaserange - L1 pseudorange */
    if (data->P[0]!=0.0&&data->L[1]!=0.0&&data->code[0]&&data->code[1]) {
        ppr=od_rtk_rtcm3e_cp_pr(data->L[1],pr1c/lam2);
        if (ppr2) *ppr2=ROUND(ppr*lam2/0.0005);
    }
    lt1=od_rtk_rtcm3e_locktime(data->time,rtcm->lltime[data->sat-1]  ,data->LLI[0]);
    lt2=od_rtk_rtcm3e_locktime(data->time,rtcm->lltime[data->sat-1]+1,data->LLI[1]);
    
    if (lock1) *lock1=od_rtk_rtcm3e_to_lock(lt1);
    if (lock2) *lock2=od_rtk_rtcm3e_to_lock(lt2);
    if (cnr1 ) *cnr1=ROUND(data->SNR[0]*SNR_UNIT/0.25);
    if (cnr2 ) *cnr2=ROUND(data->SNR[1]*SNR_UNIT/0.25);
    if (code1) *code1=od_rtk_rtcm3e_to_code1_gps(data->code[0]);
    if (code2) *code2=od_rtk_rtcm3e_to_code2_gps(data->code[1]);
}
/* generate obs field data GLONASS -------------------------------------------*/
static void od_rtk_rtcm3e_gen_obs_glo(rtcm_t *rtcm, const obsd_t *data, int fcn, int *code1,
                        int *pr1, int *ppr1, int *lock1, int *amb, int *cnr1,
                        int *code2, int *pr21, int *ppr2, int *lock2, int *cnr2)
{
    double lam1=0.0,lam2=0.0,pr1c=0.0,ppr;
    int lt1,lt2;
    
    if (fcn>=0) { /* fcn+7 */
        lam1=CLIGHT/(FREQ1_GLO+DFRQ1_GLO*(fcn-7));
        lam2=CLIGHT/(FREQ2_GLO+DFRQ2_GLO*(fcn-7));
    }
    *pr1=*amb=0;
    if (ppr1) *ppr1=0xFFF80000; /* invalid values */
    if (pr21) *pr21=0xFFFFE000;
    if (ppr2) *ppr2=0xFFF80000;
    
    /* L1 pseudorange */
    if (data->P[0]!=0.0) {
        *amb=(int)floor(data->P[0]/PRUNIT_GLO);
        *pr1=ROUND((data->P[0]-*amb*PRUNIT_GLO)/0.02);
        pr1c=*pr1*0.02+*amb*PRUNIT_GLO;
    }
    /* L1 phaserange - L1 pseudorange */
    if (data->P[0]!=0.0&&data->L[0]!=0.0&&data->code[0]&&lam1>0.0) {
        ppr=od_rtk_rtcm3e_cp_pr(data->L[0],pr1c/lam1);
        if (ppr1) *ppr1=ROUND(ppr*lam1/0.0005);
    }
    /* L2 -L1 pseudorange */
    if (data->P[0]!=0.0&&data->P[1]!=0.0&&data->code[0]&&data->code[1]&&
        fabs(data->P[1]-pr1c)<=163.82) {
        if (pr21) *pr21=ROUND((data->P[1]-pr1c)/0.02);
    }
    /* L2 phaserange - L1 pseudorange */
    if (data->P[0]!=0.0&&data->L[1]!=0.0&&data->code[0]&&data->code[1]&&
        lam2>0.0) {
        ppr=od_rtk_rtcm3e_cp_pr(data->L[1],pr1c/lam2);
        if (ppr2) *ppr2=ROUND(ppr*lam2/0.0005);
    }
    lt1=od_rtk_rtcm3e_locktime(data->time,rtcm->lltime[data->sat-1]  ,data->LLI[0]);
    lt2=od_rtk_rtcm3e_locktime(data->time,rtcm->lltime[data->sat-1]+1,data->LLI[1]);
    
    if (lock1) *lock1=od_rtk_rtcm3e_to_lock(lt1);
    if (lock2) *lock2=od_rtk_rtcm3e_to_lock(lt2);
    if (cnr1 ) *cnr1=ROUND(data->SNR[0]*SNR_UNIT/0.25);
    if (cnr2 ) *cnr2=ROUND(data->SNR[1]*SNR_UNIT/0.25);
    if (code1) *code1=od_rtk_rtcm3e_to_code1_glo(data->code[0]);
    if (code2) *code2=od_rtk_rtcm3e_to_code2_glo(data->code[1]);
}
/* encode RTCM header --------------------------------------------------------*/
static int od_rtk_rtcm3e_encode_head(int type, rtcm_t *rtcm, int sys, int sync, int nsat)
{
    double tow;
    int i=24,week,epoch;
    
    trace(4,"encode_head: type=%d sync=%d sys=%d nsat=%d\n",type,sync,sys,nsat);
    
    setbitu(rtcm->buff,i,12,type       ); i+=12; /* message no */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* ref station id */
    
    if (sys==SYS_GLO) {
        tow=time2gpst(timeadd(gpst2utc(rtcm->time),10800.0),&week);
        epoch=ROUND(fmod(tow,86400.0)/0.001);
        setbitu(rtcm->buff,i,27,epoch); i+=27; /* glonass epoch time */
    }
    else {
        tow=time2gpst(rtcm->time,&week);
        epoch=ROUND(tow/0.001);
        setbitu(rtcm->buff,i,30,epoch); i+=30; /* gps epoch time */
    }
    setbitu(rtcm->buff,i, 1,sync); i+= 1; /* synchronous gnss flag */
    setbitu(rtcm->buff,i, 5,nsat); i+= 5; /* no of satellites */
    setbitu(rtcm->buff,i, 1,0   ); i+= 1; /* smoothing indicator */
    setbitu(rtcm->buff,i, 3,0   ); i+= 3; /* smoothing interval */
    return i;
}
/* encode type 1001: basic L1-only GPS RTK observables -----------------------*/
static int od_rtk_rtcm3e_encode_type1001(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sys,prn;
    int code1,pr1,ppr1,lock1,amb;
    
    trace(3,"encode_type1001: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1001,rtcm,SYS_GPS,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        
        if (sys==SYS_SBS) prn-=80; /* 40-58: sbas 120-138 */
        
        /* generate obs field data gps */
        od_rtk_rtcm3e_gen_obs_gps(rtcm,rtcm->obs.data+j,&code1,&pr1,&ppr1,&lock1,&amb,NULL,
                    NULL,NULL,NULL,NULL,NULL);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i,24,pr1  ); i+=24;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1002: extended L1-only GPS RTK observables --------------------*/
static int od_rtk_rtcm3e_encode_type1002(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sys,prn;
    int code1,pr1,ppr1,lock1,amb,cnr1;
    
    trace(3,"encode_type1002: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1002,rtcm,SYS_GPS,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        
        if (sys==SYS_SBS) prn-=80; /* 40-58: sbas 120-138 */
        
        /* generate obs field data gps */
        od_rtk_rtcm3e_gen_obs_gps(rtcm,rtcm->obs.data+j,&code1,&pr1,&ppr1,&lock1,&amb,&cnr1,
                    NULL,NULL,NULL,NULL,NULL);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i,24,pr1  ); i+=24;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
        setbitu(rtcm->buff,i, 8,amb  ); i+= 8;
        setbitu(rtcm->buff,i, 8,cnr1 ); i+= 8;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1003: basic L1&L2 GPS RTK observables -------------------------*/
static int od_rtk_rtcm3e_encode_type1003(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sys,prn;
    int code1,pr1,ppr1,lock1,amb,code2,pr21,ppr2,lock2;
    
    trace(3,"encode_type1003: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1003,rtcm,SYS_GPS,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        
        if (sys==SYS_SBS) prn-=80; /* 40-58: sbas 120-138 */
        
        /* generate obs field data gps */
        od_rtk_rtcm3e_gen_obs_gps(rtcm,rtcm->obs.data+j,&code1,&pr1,&ppr1,&lock1,&amb,
                    NULL,&code2,&pr21,&ppr2,&lock2,NULL);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i,24,pr1  ); i+=24;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
        setbitu(rtcm->buff,i, 2,code2); i+= 2;
        setbits(rtcm->buff,i,14,pr21 ); i+=14;
        setbits(rtcm->buff,i,20,ppr2 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock2); i+= 7;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1004: extended L1&L2 GPS RTK observables ----------------------*/
static int od_rtk_rtcm3e_encode_type1004(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sys,prn;
    int code1,pr1,ppr1,lock1,amb,cnr1,code2,pr21,ppr2,lock2,cnr2;
    
    trace(3,"encode_type1004: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1004,rtcm,SYS_GPS,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sys=satsys(rtcm->obs.data[j].sat,&prn);
        if (!(sys&(SYS_GPS|SYS_SBS))) continue;
        
        if (sys==SYS_SBS) prn-=80; /* 40-58: sbas 120-138 */
        
        /* generate obs field data gps */
        od_rtk_rtcm3e_gen_obs_gps(rtcm,rtcm->obs.data+j,&code1,&pr1,&ppr1,&lock1,&amb,
                    &cnr1,&code2,&pr21,&ppr2,&lock2,&cnr2);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i,24,pr1  ); i+=24;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
        setbitu(rtcm->buff,i, 8,amb  ); i+= 8;
        setbitu(rtcm->buff,i, 8,cnr1 ); i+= 8;
        setbitu(rtcm->buff,i, 2,code2); i+= 2;
        setbits(rtcm->buff,i,14,pr21 ); i+=14;
        setbits(rtcm->buff,i,20,ppr2 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock2); i+= 7;
        setbitu(rtcm->buff,i, 8,cnr2 ); i+= 8;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1005: stationary RTK reference station ARP --------------------*/
static int od_rtk_rtcm3e_encode_type1005(rtcm_t *rtcm, int sync)
{
    double *p=rtcm->sta.pos;
    int i=24;
    
    trace(3,"encode_type1005: sync=%d\n",sync);
    
    setbitu(rtcm->buff,i,12,1005       ); i+=12; /* message no */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* ref station id */
    setbitu(rtcm->buff,i, 6,0          ); i+= 6; /* itrf realization year */
    setbitu(rtcm->buff,i, 1,1          ); i+= 1; /* gps indicator */
    setbitu(rtcm->buff,i, 1,1          ); i+= 1; /* glonass indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* galileo indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* ref station indicator */
    od_rtk_rtcm3e_set38bits(rtcm->buff,i,p[0]/0.0001 ); i+=38; /* antenna ref point ecef-x */
    setbitu(rtcm->buff,i, 1,1          ); i+= 1; /* oscillator indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* reserved */
    od_rtk_rtcm3e_set38bits(rtcm->buff,i,p[1]/0.0001 ); i+=38; /* antenna ref point ecef-y */
    setbitu(rtcm->buff,i, 2,0          ); i+= 2; /* quarter cycle indicator */
    od_rtk_rtcm3e_set38bits(rtcm->buff,i,p[2]/0.0001 ); i+=38; /* antenna ref point ecef-z */
    rtcm->nbit=i;
    return 1;
}
/* encode type 1006: stationary RTK reference station ARP with height --------*/
static int od_rtk_rtcm3e_encode_type1006(rtcm_t *rtcm, int sync)
{
    double *p=rtcm->sta.pos;
    int i=24,hgt=0;
    
    trace(3,"encode_type1006: sync=%d\n",sync);
    
    if (0.0<=rtcm->sta.hgt&&rtcm->sta.hgt<=6.5535) {
        hgt=ROUND(rtcm->sta.hgt/0.0001);
    }
    else {
        trace(2,"antenna height error: h=%.4f\n",rtcm->sta.hgt);
    }
    setbitu(rtcm->buff,i,12,1006       ); i+=12; /* message no */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* ref station id */
    setbitu(rtcm->buff,i, 6,0          ); i+= 6; /* itrf realization year */
    setbitu(rtcm->buff,i, 1,1          ); i+= 1; /* gps indicator */
    setbitu(rtcm->buff,i, 1,1          ); i+= 1; /* glonass indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* galileo indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* ref station indicator */
    od_rtk_rtcm3e_set38bits(rtcm->buff,i,p[0]/0.0001 ); i+=38; /* antenna ref point ecef-x */
    setbitu(rtcm->buff,i, 1,1          ); i+= 1; /* oscillator indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* reserved */
    od_rtk_rtcm3e_set38bits(rtcm->buff,i,p[1]/0.0001 ); i+=38; /* antenna ref point ecef-y */
    setbitu(rtcm->buff,i, 2,0          ); i+= 2; /* quarter cycle indicator */
    od_rtk_rtcm3e_set38bits(rtcm->buff,i,p[2]/0.0001 ); i+=38; /* antenna ref point ecef-z */
    setbitu(rtcm->buff,i,16,hgt        ); i+=16; /* antenna height */
    rtcm->nbit=i;
    return 1;
}
/* encode type 1007: antenna descriptor --------------------------------------*/
static int od_rtk_rtcm3e_encode_type1007(rtcm_t *rtcm, int sync)
{
    int i=24,j,antsetup=rtcm->sta.antsetup;
    int n=MIN(strlen(rtcm->sta.antdes),31);
    
    trace(3,"encode_type1007: sync=%d\n",sync);
    
    setbitu(rtcm->buff,i,12,1007       ); i+=12; /* message no */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* ref station id */
    
    /* antenna descriptor */
    setbitu(rtcm->buff,i,8,n); i+=8;
    for (j=0;j<n;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.antdes[j]); i+=8;
    }
    setbitu(rtcm->buff,i,8,antsetup); i+=8; /* antetnna setup id */
    rtcm->nbit=i;
    return 1;
}
/* encode type 1008: antenna descriptor & serial number ----------------------*/
static int od_rtk_rtcm3e_encode_type1008(rtcm_t *rtcm, int sync)
{
    int i=24,j,antsetup=rtcm->sta.antsetup;
    int n=MIN(strlen(rtcm->sta.antdes),31);
    int m=MIN(strlen(rtcm->sta.antsno),31);
    
    trace(3,"encode_type1008: sync=%d\n",sync);
    
    setbitu(rtcm->buff,i,12,1008       ); i+=12; /* message no */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* ref station id */
    
    /* antenna descriptor */
    setbitu(rtcm->buff,i,8,n); i+=8;
    for (j=0;j<n;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.antdes[j]); i+=8;
    }
    setbitu(rtcm->buff,i,8,antsetup); i+=8; /* antenna setup id */
    
    /* antenna serial number */
    setbitu(rtcm->buff,i,8,m); i+=8;
    for (j=0;j<m;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.antsno[j]); i+=8;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1009: basic L1-only GLONASS RTK observables -------------------*/
static int od_rtk_rtcm3e_encode_type1009(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sat,prn,fcn;
    int code1,pr1,ppr1,lock1,amb;
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1009,rtcm,SYS_GLO,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        
        /* generate obs field data glonass */
        od_rtk_rtcm3e_gen_obs_glo(rtcm,rtcm->obs.data+j,fcn,&code1,&pr1,&ppr1,&lock1,&amb,
                    NULL,NULL,NULL,NULL,NULL,NULL);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i, 5,fcn  ); i+= 5; /* fcn+7 */
        setbitu(rtcm->buff,i,25,pr1  ); i+=25;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1010: extended L1-only GLONASS RTK observables ----------------*/
static int od_rtk_rtcm3e_encode_type1010(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sat,prn,fcn;
    int code1,pr1,ppr1,lock1,amb,cnr1;
    
    trace(3,"encode_type1010: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1010,rtcm,SYS_GLO,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        
        /* generate obs field data glonass */
        od_rtk_rtcm3e_gen_obs_glo(rtcm,rtcm->obs.data+j,fcn,&code1,&pr1,&ppr1,&lock1,&amb,
                    &cnr1,NULL,NULL,NULL,NULL,NULL);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i, 5,fcn  ); i+= 5; /* fcn+7 */
        setbitu(rtcm->buff,i,25,pr1  ); i+=25;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
        setbitu(rtcm->buff,i, 7,amb  ); i+= 7;
        setbitu(rtcm->buff,i, 8,cnr1 ); i+= 8;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1011: basic  L1&L2 GLONASS RTK observables --------------------*/
static int od_rtk_rtcm3e_encode_type1011(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sat,prn,fcn;
    int code1,pr1,ppr1,lock1,amb,code2,pr21,ppr2,lock2;
    
    trace(3,"encode_type1011: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1011,rtcm,SYS_GLO,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        
        /* generate obs field data glonass */
        od_rtk_rtcm3e_gen_obs_glo(rtcm,rtcm->obs.data+j,fcn,&code1,&pr1,&ppr1,&lock1,&amb,
                    NULL,&code2,&pr21,&ppr2,&lock2,NULL);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i, 5,fcn  ); i+= 5; /* fcn+7 */
        setbitu(rtcm->buff,i,25,pr1  ); i+=25;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
        setbitu(rtcm->buff,i, 2,code2); i+= 2;
        setbits(rtcm->buff,i,14,pr21 ); i+=14;
        setbits(rtcm->buff,i,20,ppr2 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock2); i+= 7;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1012: extended L1&L2 GLONASS RTK observables ------------------*/
static int od_rtk_rtcm3e_encode_type1012(rtcm_t *rtcm, int sync)
{
    int i,j,nsat=0,sat,prn,fcn;
    int code1,pr1,ppr1,lock1,amb,cnr1,code2,pr21,ppr2,lock2,cnr2;
    
    trace(3,"encode_type1012: sync=%d\n",sync);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue;  /* fcn+7 */
        nsat++;
    }
    /* encode header */
    i=od_rtk_rtcm3e_encode_head(1012,rtcm,SYS_GLO,sync,nsat);
    
    for (j=0;j<rtcm->obs.n&&nsat<MAXOBS;j++) {
        sat=rtcm->obs.data[j].sat;
        if (satsys(sat,&prn)!=SYS_GLO) continue;
        if ((fcn=od_rtk_rtcm3e_fcn_glo(sat,rtcm))<0) continue; /* fcn+7 */
        
        /* generate obs field data glonass */
        od_rtk_rtcm3e_gen_obs_glo(rtcm,rtcm->obs.data+j,fcn,&code1,&pr1,&ppr1,&lock1,&amb,
                    &cnr1,&code2,&pr21,&ppr2,&lock2,&cnr2);
        
        setbitu(rtcm->buff,i, 6,prn  ); i+= 6;
        setbitu(rtcm->buff,i, 1,code1); i+= 1;
        setbitu(rtcm->buff,i, 5,fcn  ); i+= 5; /* fcn+7 */
        setbitu(rtcm->buff,i,25,pr1  ); i+=25;
        setbits(rtcm->buff,i,20,ppr1 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock1); i+= 7;
        setbitu(rtcm->buff,i, 7,amb  ); i+= 7;
        setbitu(rtcm->buff,i, 8,cnr1 ); i+= 8;
        setbitu(rtcm->buff,i, 2,code2); i+= 2;
        setbits(rtcm->buff,i,14,pr21 ); i+=14;
        setbits(rtcm->buff,i,20,ppr2 ); i+=20;
        setbitu(rtcm->buff,i, 7,lock2); i+= 7;
        setbitu(rtcm->buff,i, 8,cnr2 ); i+= 8;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1019: GPS ephemerides -----------------------------------------*/
static int od_rtk_rtcm3e_encode_type1019(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,tgd;
    
    trace(3,"encode_type1019: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_GPS) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1;
    if (eph->sat!=rtcm->ephsat) return 0;
    week=eph->week%1024;
    toe  =ROUND(eph->toes/16.0);
    toc  =ROUND(time2gpst(eph->toc,NULL)/16.0);
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_43/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_43/SC2RAD);
    crs  =ROUND(eph->crs/P2_5 );
    crc  =ROUND(eph->crc/P2_5 );
    cus  =ROUND(eph->cus/P2_29);
    cuc  =ROUND(eph->cuc/P2_29);
    cis  =ROUND(eph->cis/P2_29);
    cic  =ROUND(eph->cic/P2_29);
    af0  =ROUND(eph->f0 /P2_31);
    af1  =ROUND(eph->f1 /P2_43);
    af2  =ROUND(eph->f2 /P2_55);
    tgd  =ROUND(eph->tgd[0]/P2_31);
    
    setbitu(rtcm->buff,i,12,1019     ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i,10,week     ); i+=10;
    setbitu(rtcm->buff,i, 4,eph->sva ); i+= 4;
    setbitu(rtcm->buff,i, 2,eph->code); i+= 2;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbitu(rtcm->buff,i, 8,eph->iode); i+= 8;
    setbitu(rtcm->buff,i,16,toc      ); i+=16;
    setbits(rtcm->buff,i, 8,af2      ); i+= 8;
    setbits(rtcm->buff,i,16,af1      ); i+=16;
    setbits(rtcm->buff,i,22,af0      ); i+=22;
    setbitu(rtcm->buff,i,10,eph->iodc); i+=10;
    setbits(rtcm->buff,i,16,crs      ); i+=16;
    setbits(rtcm->buff,i,16,deln     ); i+=16;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbits(rtcm->buff,i,16,cuc      ); i+=16;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbits(rtcm->buff,i,16,cus      ); i+=16;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbitu(rtcm->buff,i,16,toe      ); i+=16;
    setbits(rtcm->buff,i,16,cic      ); i+=16;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,16,cis      ); i+=16;
    setbits(rtcm->buff,i,32,i0       ); i+=32;
    setbits(rtcm->buff,i,16,crc      ); i+=16;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,24,OMGd     ); i+=24;
    setbits(rtcm->buff,i, 8,tgd      ); i+= 8;
    setbitu(rtcm->buff,i, 6,eph->svh ); i+= 6;
    setbitu(rtcm->buff,i, 1,eph->flag); i+= 1;
    setbitu(rtcm->buff,i, 1,eph->fit>0.0?0:1); i+=1;
    rtcm->nbit=i;
    return 1;
}
/* encode type 1020: GLONASS ephemerides -------------------------------------*/
static int od_rtk_rtcm3e_encode_type1020(rtcm_t *rtcm, int sync)
{
    geph_t *geph;
    gtime_t time;
    double ep[6];
    int i=24,j,prn,tk_h,tk_m,tk_s,tb,pos[3],vel[3],acc[3],gamn,taun,dtaun;
    int fcn,NT;
    
    trace(3,"encode_type1020: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_GLO) return 0;
    geph=rtcm->nav.geph+prn-1;
    if (geph->sat!=rtcm->ephsat) return 0;
    fcn=geph->frq+7;
    
    /* time of frame within day (utc(su) + 3 hr) */
    time=timeadd(gpst2utc(geph->tof),10800.0);
    time2epoch(time,ep);
    tk_h=(int)ep[3];
    tk_m=(int)ep[4];
    tk_s=ROUND(ep[5]/30.0);
    
    /* # of days since jan 1 in leap year */
    ep[0]=floor(ep[0]/4.0)*4.0; ep[1]=ep[2]=1.0;
    ep[3]=ep[4]=ep[5]=0.0;
    NT=(int)floor(timediff(time,epoch2time(ep))/86400.+1.0);
    
    /* index of time interval within day (utc(su) + 3 hr) */
    time=timeadd(gpst2utc(geph->toe),10800.0);
    time2epoch(time,ep);
    tb=ROUND((ep[3]*3600.0+ep[4]*60.0+ep[5])/900.0);
    
    for (j=0;j<3;j++) {
        pos[j]=ROUND(geph->pos[j]/P2_11/1E3);
        vel[j]=ROUND(geph->vel[j]/P2_20/1E3);
        acc[j]=ROUND(geph->acc[j]/P2_30/1E3);
    }
    gamn =ROUND(geph->gamn /P2_40);
    taun =ROUND(geph->taun /P2_30);
    dtaun=ROUND(geph->dtaun/P2_30);
    
    setbitu(rtcm->buff,i,12,1020     ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i, 5,fcn      ); i+= 5;
    setbitu(rtcm->buff,i, 4,0        ); i+= 4; /* almanac health,P1 */
    setbitu(rtcm->buff,i, 5,tk_h     ); i+= 5;
    setbitu(rtcm->buff,i, 6,tk_m     ); i+= 6;
    setbitu(rtcm->buff,i, 1,tk_s     ); i+= 1;
    setbitu(rtcm->buff,i, 1,geph->svh); i+= 1; /* Bn */
    setbitu(rtcm->buff,i, 1,0        ); i+= 1; /* P2 */
    setbitu(rtcm->buff,i, 7,tb       ); i+= 7;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,24,vel[0]   ); i+=24;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,27,pos[0]   ); i+=27;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i, 5,acc[0]   ); i+= 5;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,24,vel[1]   ); i+=24;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,27,pos[1]   ); i+=27;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i, 5,acc[1]   ); i+= 5;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,24,vel[2]   ); i+=24;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,27,pos[2]   ); i+=27;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i, 5,acc[2]   ); i+= 5;
    setbitu(rtcm->buff,i, 1,0        ); i+= 1; /* P3 */
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,11,gamn     ); i+=11;
    setbitu(rtcm->buff,i, 3,0        ); i+= 3; /* P,ln */
    od_rtk_rtcm3e_setbitg(rtcm->buff,i,22,taun     ); i+=22;
    od_rtk_rtcm3e_setbitg(rtcm->buff,i, 5,dtaun    ); i+= 5;
    setbitu(rtcm->buff,i, 5,geph->age); i+= 5; /* En */
    setbitu(rtcm->buff,i, 1,0        ); i+= 1; /* P4 */
    setbitu(rtcm->buff,i, 4,0        ); i+= 4; /* FT */
    setbitu(rtcm->buff,i,11,NT       ); i+=11;
    setbitu(rtcm->buff,i, 2,0        ); i+= 2; /* M */
    setbitu(rtcm->buff,i, 1,0        ); i+= 1; /* flag for additional data */
    setbitu(rtcm->buff,i,11,0        ); i+=11; /* NA */
    setbitu(rtcm->buff,i,32,0        ); i+=32; /* tauc */
    setbitu(rtcm->buff,i, 5,0        ); i+= 5; /* N4 */
    setbitu(rtcm->buff,i,22,0        ); i+=22; /* taugps */
    setbitu(rtcm->buff,i, 1,0        ); i+= 1; /* ln */
    setbitu(rtcm->buff,i, 7,0        ); i+= 7;
    rtcm->nbit=i;
    return 1;
}
/* encode type 1033: receiver and antenna descriptor -------------------------*/
static int od_rtk_rtcm3e_encode_type1033(rtcm_t *rtcm, int sync)
{
    int i=24,j,antsetup=rtcm->sta.antsetup;
    int n=MIN(strlen(rtcm->sta.antdes ),31);
    int m=MIN(strlen(rtcm->sta.antsno ),31);
    int I=MIN(strlen(rtcm->sta.rectype),31);
    int J=MIN(strlen(rtcm->sta.recver ),31);
    int K=MIN(strlen(rtcm->sta.recsno ),31);
    
    trace(3,"encode_type1033: sync=%d\n",sync);
    
    setbitu(rtcm->buff,i,12,1033       ); i+=12;
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12;
    
    setbitu(rtcm->buff,i,8,n); i+= 8;
    for (j=0;j<n;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.antdes[j]); i+=8;
    }
    setbitu(rtcm->buff,i,8,antsetup); i+= 8;
    
    setbitu(rtcm->buff,i,8,m); i+= 8;
    for (j=0;j<m;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.antsno[j]); i+=8;
    }
    setbitu(rtcm->buff,i,8,I); i+= 8;
    for (j=0;j<I;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.rectype[j]); i+=8;
    }
    setbitu(rtcm->buff,i,8,J); i+= 8;
    for (j=0;j<J;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.recver[j]); i+=8;
    }
    setbitu(rtcm->buff,i,8,K); i+= 8;
    for (j=0;j<K;j++) {
        setbitu(rtcm->buff,i,8,rtcm->sta.recsno[j]); i+=8;
    }
    rtcm->nbit=i;
    return 1;
}
/* encode type 1041: NavIC/IRNSS ephemerides ---------------------------------*/
static int od_rtk_rtcm3e_encode_type1041(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,tgd;
    
    trace(3,"encode_type1041: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_IRN) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1;
    if (eph->sat!=rtcm->ephsat) return 0;
    week=eph->week%1024;
    toe  =ROUND(eph->toes/16.0);
    toc  =ROUND(time2gpst(eph->toc,NULL)/16.0);
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_41/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_41/SC2RAD);
    crs  =ROUND(eph->crs/0.0625);
    crc  =ROUND(eph->crc/0.0625);
    cus  =ROUND(eph->cus/P2_28);
    cuc  =ROUND(eph->cuc/P2_28);
    cis  =ROUND(eph->cis/P2_28);
    cic  =ROUND(eph->cic/P2_28);
    af0  =ROUND(eph->f0 /P2_31);
    af1  =ROUND(eph->f1 /P2_43);
    af2  =ROUND(eph->f2 /P2_55);
    tgd  =ROUND(eph->tgd[0]/P2_31);
    
    setbitu(rtcm->buff,i,12,1041     ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i,10,week     ); i+=10;
    setbits(rtcm->buff,i,22,af0      ); i+=22;
    setbits(rtcm->buff,i,16,af1      ); i+=16;
    setbits(rtcm->buff,i, 8,af2      ); i+= 8;
    setbitu(rtcm->buff,i, 4,eph->sva ); i+= 4;
    setbitu(rtcm->buff,i,16,toc      ); i+=16;
    setbits(rtcm->buff,i, 8,tgd      ); i+= 8;
    setbits(rtcm->buff,i,22,deln     ); i+=22;
    setbitu(rtcm->buff,i, 8,eph->iode); i+= 8+10; /* IODEC */
    setbitu(rtcm->buff,i, 2,eph->svh ); i+= 2; /* L5+Sflag */
    setbits(rtcm->buff,i,15,cuc      ); i+=15;
    setbits(rtcm->buff,i,15,cus      ); i+=15;
    setbits(rtcm->buff,i,15,cic      ); i+=15;
    setbits(rtcm->buff,i,15,cis      ); i+=15;
    setbits(rtcm->buff,i,15,crc      ); i+=15;
    setbits(rtcm->buff,i,15,crs      ); i+=15;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbitu(rtcm->buff,i,16,toe      ); i+=16;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,22,OMGd     ); i+=22;
    setbits(rtcm->buff,i,32,i0       ); i+=32+4;
    rtcm->nbit=i;
    return 1;
}
/* encode type 1044: QZSS ephemerides ----------------------------------------*/
static int od_rtk_rtcm3e_encode_type1044(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,tgd;
    
    trace(3,"encode_type1044: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_QZS) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1;
    if (eph->sat!=rtcm->ephsat) return 0;
    week=eph->week%1024;
    toe  =ROUND(eph->toes/16.0);
    toc  =ROUND(time2gpst(eph->toc,NULL)/16.0);
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_43/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_43/SC2RAD);
    crs  =ROUND(eph->crs/P2_5 );
    crc  =ROUND(eph->crc/P2_5 );
    cus  =ROUND(eph->cus/P2_29);
    cuc  =ROUND(eph->cuc/P2_29);
    cis  =ROUND(eph->cis/P2_29);
    cic  =ROUND(eph->cic/P2_29);
    af0  =ROUND(eph->f0 /P2_31);
    af1  =ROUND(eph->f1 /P2_43);
    af2  =ROUND(eph->f2 /P2_55);
    tgd  =ROUND(eph->tgd[0]/P2_31);
    
    setbitu(rtcm->buff,i,12,1044     ); i+=12;
    setbitu(rtcm->buff,i, 4,prn-192  ); i+= 4;
    setbitu(rtcm->buff,i,16,toc      ); i+=16;
    setbits(rtcm->buff,i, 8,af2      ); i+= 8;
    setbits(rtcm->buff,i,16,af1      ); i+=16;
    setbits(rtcm->buff,i,22,af0      ); i+=22;
    setbitu(rtcm->buff,i, 8,eph->iode); i+= 8;
    setbits(rtcm->buff,i,16,crs      ); i+=16;
    setbits(rtcm->buff,i,16,deln     ); i+=16;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbits(rtcm->buff,i,16,cuc      ); i+=16;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbits(rtcm->buff,i,16,cus      ); i+=16;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbitu(rtcm->buff,i,16,toe      ); i+=16;
    setbits(rtcm->buff,i,16,cic      ); i+=16;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,16,cis      ); i+=16;
    setbits(rtcm->buff,i,32,i0       ); i+=32;
    setbits(rtcm->buff,i,16,crc      ); i+=16;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,24,OMGd     ); i+=24;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbitu(rtcm->buff,i, 2,eph->code); i+= 2;
    setbitu(rtcm->buff,i,10,week     ); i+=10;
    setbitu(rtcm->buff,i, 4,eph->sva ); i+= 4;
    setbitu(rtcm->buff,i, 6,eph->svh ); i+= 6;
    setbits(rtcm->buff,i, 8,tgd      ); i+= 8;
    setbitu(rtcm->buff,i,10,eph->iodc); i+=10;
    setbitu(rtcm->buff,i, 1,eph->fit==2.0?0:1); i+=1;
    rtcm->nbit=i;
    return 1;
}
/* encode type 1045: Galileo F/NAV satellite ephemerides ---------------------*/
static int od_rtk_rtcm3e_encode_type1045(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,bgd1,bgd2,oshs,osdvs;
    
    trace(3,"encode_type1045: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_GAL) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1+MAXSAT; /* F/NAV */
    if (eph->sat!=rtcm->ephsat) return 0;
    week=(eph->week-1024)%4096; /* gst-week = gal-week - 1024 */
    toe  =ROUND(eph->toes/60.0);
    toc  =ROUND(time2gpst(eph->toc,NULL)/60.0);
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_43/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_43/SC2RAD);
    crs  =ROUND(eph->crs/P2_5 );
    crc  =ROUND(eph->crc/P2_5 );
    cus  =ROUND(eph->cus/P2_29);
    cuc  =ROUND(eph->cuc/P2_29);
    cis  =ROUND(eph->cis/P2_29);
    cic  =ROUND(eph->cic/P2_29);
    af0  =ROUND(eph->f0 /P2_34);
    af1  =ROUND(eph->f1 /P2_46);
    af2  =ROUND(eph->f2 /P2_59);
    bgd1 =ROUND(eph->tgd[0]/P2_32); /* E5a/E1 */
    bgd2 =ROUND(eph->tgd[1]/P2_32); /* E5b/E1 */
    oshs =(eph->svh>>4)&3;          /* E5a SVH */
    osdvs=(eph->svh>>3)&1;          /* E5a DVS */
    setbitu(rtcm->buff,i,12,1045     ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i,12,week     ); i+=12;
    setbitu(rtcm->buff,i,10,eph->iode); i+=10;
    setbitu(rtcm->buff,i, 8,eph->sva ); i+= 8;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbitu(rtcm->buff,i,14,toc      ); i+=14;
    setbits(rtcm->buff,i, 6,af2      ); i+= 6;
    setbits(rtcm->buff,i,21,af1      ); i+=21;
    setbits(rtcm->buff,i,31,af0      ); i+=31;
    setbits(rtcm->buff,i,16,crs      ); i+=16;
    setbits(rtcm->buff,i,16,deln     ); i+=16;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbits(rtcm->buff,i,16,cuc      ); i+=16;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbits(rtcm->buff,i,16,cus      ); i+=16;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbitu(rtcm->buff,i,14,toe      ); i+=14;
    setbits(rtcm->buff,i,16,cic      ); i+=16;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,16,cis      ); i+=16;
    setbits(rtcm->buff,i,32,i0       ); i+=32;
    setbits(rtcm->buff,i,16,crc      ); i+=16;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,24,OMGd     ); i+=24;
    setbits(rtcm->buff,i,10,bgd1     ); i+=10;
    setbitu(rtcm->buff,i, 2,oshs     ); i+= 2; /* E5a SVH */
    setbitu(rtcm->buff,i, 1,osdvs    ); i+= 1; /* E5a DVS */
    setbitu(rtcm->buff,i, 7,0        ); i+= 7; /* reserved */
    rtcm->nbit=i;
    return 1;
}
/* encode type 1046: Galileo I/NAV satellite ephemerides ---------------------*/
static int od_rtk_rtcm3e_encode_type1046(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,bgd1,bgd2,oshs1,osdvs1,oshs2,osdvs2;
    
    trace(3,"encode_type1046: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_GAL) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1; /* I/NAV */
    if (eph->sat!=rtcm->ephsat) return 0;
    week=(eph->week-1024)%4096; /* gst-week = gal-week - 1024 */
    toe  =ROUND(eph->toes/60.0);
    toc  =ROUND(time2gpst(eph->toc,NULL)/60.0);
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_43/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_43/SC2RAD);
    crs  =ROUND(eph->crs/P2_5 );
    crc  =ROUND(eph->crc/P2_5 );
    cus  =ROUND(eph->cus/P2_29);
    cuc  =ROUND(eph->cuc/P2_29);
    cis  =ROUND(eph->cis/P2_29);
    cic  =ROUND(eph->cic/P2_29);
    af0  =ROUND(eph->f0 /P2_34);
    af1  =ROUND(eph->f1 /P2_46);
    af2  =ROUND(eph->f2 /P2_59);
    bgd1 =ROUND(eph->tgd[0]/P2_32); /* E5a/E1 */
    bgd2 =ROUND(eph->tgd[1]/P2_32); /* E5b/E1 */
    oshs1 =(eph->svh>>7)&3;         /* E5b SVH */
    osdvs1=(eph->svh>>6)&1;         /* E5b DVS */
    oshs2 =(eph->svh>>1)&3;         /* E1 SVH */
    osdvs2=(eph->svh>>0)&1;         /* E1 DVS */
    setbitu(rtcm->buff,i,12,1046     ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i,12,week     ); i+=12;
    setbitu(rtcm->buff,i,10,eph->iode); i+=10;
    setbitu(rtcm->buff,i, 8,eph->sva ); i+= 8;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbitu(rtcm->buff,i,14,toc      ); i+=14;
    setbits(rtcm->buff,i, 6,af2      ); i+= 6;
    setbits(rtcm->buff,i,21,af1      ); i+=21;
    setbits(rtcm->buff,i,31,af0      ); i+=31;
    setbits(rtcm->buff,i,16,crs      ); i+=16;
    setbits(rtcm->buff,i,16,deln     ); i+=16;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbits(rtcm->buff,i,16,cuc      ); i+=16;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbits(rtcm->buff,i,16,cus      ); i+=16;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbitu(rtcm->buff,i,14,toe      ); i+=14;
    setbits(rtcm->buff,i,16,cic      ); i+=16;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,16,cis      ); i+=16;
    setbits(rtcm->buff,i,32,i0       ); i+=32;
    setbits(rtcm->buff,i,16,crc      ); i+=16;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,24,OMGd     ); i+=24;
    setbits(rtcm->buff,i,10,bgd1     ); i+=10;
    setbits(rtcm->buff,i,10,bgd2     ); i+=10;
    setbitu(rtcm->buff,i, 2,oshs1    ); i+= 2; /* E5b SVH */
    setbitu(rtcm->buff,i, 1,osdvs1   ); i+= 1; /* E5b DVS */
    setbitu(rtcm->buff,i, 2,oshs2    ); i+= 2; /* E1 SVH */
    setbitu(rtcm->buff,i, 1,osdvs2   ); i+= 1; /* E1 DVS */
    rtcm->nbit=i;
    return 1;
}
/* encode type 1042: Beidou ephemerides --------------------------------------*/
static int od_rtk_rtcm3e_encode_type1042(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,tgd1,tgd2;
    
    trace(3,"encode_type1042: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_CMP) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1;
    if (eph->sat!=rtcm->ephsat) return 0;
    week =eph->week%8192;
    toe  =ROUND(eph->toes/8.0);
    toc  =ROUND(time2bdt(gpst2bdt(eph->toc),NULL)/8.0); /* gpst -> bdt */
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_43/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_43/SC2RAD);
    crs  =ROUND(eph->crs/P2_6 );
    crc  =ROUND(eph->crc/P2_6 );
    cus  =ROUND(eph->cus/P2_31);
    cuc  =ROUND(eph->cuc/P2_31);
    cis  =ROUND(eph->cis/P2_31);
    cic  =ROUND(eph->cic/P2_31);
    af0  =ROUND(eph->f0 /P2_33);
    af1  =ROUND(eph->f1 /P2_50);
    af2  =ROUND(eph->f2 /P2_66);
    tgd1 =ROUND(eph->tgd[0]/1E-10);
    tgd2 =ROUND(eph->tgd[1]/1E-10);
    
    setbitu(rtcm->buff,i,12,1042     ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i,13,week     ); i+=13;
    setbitu(rtcm->buff,i, 4,eph->sva ); i+= 4;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbitu(rtcm->buff,i, 5,eph->iode); i+= 5;
    setbitu(rtcm->buff,i,17,toc      ); i+=17;
    setbits(rtcm->buff,i,11,af2      ); i+=11;
    setbits(rtcm->buff,i,22,af1      ); i+=22;
    setbits(rtcm->buff,i,24,af0      ); i+=24;
    setbitu(rtcm->buff,i, 5,eph->iodc); i+= 5;
    setbits(rtcm->buff,i,18,crs      ); i+=18;
    setbits(rtcm->buff,i,16,deln     ); i+=16;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbits(rtcm->buff,i,18,cuc      ); i+=18;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbits(rtcm->buff,i,18,cus      ); i+=18;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbitu(rtcm->buff,i,17,toe      ); i+=17;
    setbits(rtcm->buff,i,18,cic      ); i+=18;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,18,cis      ); i+=18;
    setbits(rtcm->buff,i,32,i0       ); i+=32;
    setbits(rtcm->buff,i,18,crc      ); i+=18;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,24,OMGd     ); i+=24;
    setbits(rtcm->buff,i,10,tgd1     ); i+=10;
    setbits(rtcm->buff,i,10,tgd2     ); i+=10;
    setbitu(rtcm->buff,i, 1,eph->svh ); i+= 1;
    rtcm->nbit=i;
    return 1;
}
/* encode type 63: Beidou ephemerides (RTCM draft) ---------------------------*/
static int od_rtk_rtcm3e_encode_type63(rtcm_t *rtcm, int sync)
{
    eph_t *eph;
    uint32_t sqrtA,e;
    int i=24,prn,week,toe,toc,i0,OMG0,omg,M0,deln,idot,OMGd,crs,crc;
    int cus,cuc,cis,cic,af0,af1,af2,tgd1,tgd2;
    
    trace(3,"encode_type63: sync=%d\n",sync);
    
    if (satsys(rtcm->ephsat,&prn)!=SYS_CMP) return 0;
    eph=rtcm->nav.eph+rtcm->ephsat-1;
    if (eph->sat!=rtcm->ephsat) return 0;
    week =eph->week%8192;
    toe  =ROUND(eph->toes/8.0);
    toc  =ROUND(time2bdt(gpst2bdt(eph->toc),NULL)/8.0); /* gpst -> bdt */
    sqrtA=ROUND_U(sqrt(eph->A)/P2_19);
    e    =ROUND_U(eph->e/P2_33);
    i0   =ROUND(eph->i0  /P2_31/SC2RAD);
    OMG0 =ROUND(eph->OMG0/P2_31/SC2RAD);
    omg  =ROUND(eph->omg /P2_31/SC2RAD);
    M0   =ROUND(eph->M0  /P2_31/SC2RAD);
    deln =ROUND(eph->deln/P2_43/SC2RAD);
    idot =ROUND(eph->idot/P2_43/SC2RAD);
    OMGd =ROUND(eph->OMGd/P2_43/SC2RAD);
    crs  =ROUND(eph->crs/P2_6 );
    crc  =ROUND(eph->crc/P2_6 );
    cus  =ROUND(eph->cus/P2_31);
    cuc  =ROUND(eph->cuc/P2_31);
    cis  =ROUND(eph->cis/P2_31);
    cic  =ROUND(eph->cic/P2_31);
    af0  =ROUND(eph->f0 /P2_33);
    af1  =ROUND(eph->f1 /P2_50);
    af2  =ROUND(eph->f2 /P2_66);
    tgd1 =ROUND(eph->tgd[0]/1E-10);
    tgd2 =ROUND(eph->tgd[1]/1E-10);
    
    setbitu(rtcm->buff,i,12,63       ); i+=12;
    setbitu(rtcm->buff,i, 6,prn      ); i+= 6;
    setbitu(rtcm->buff,i,13,week     ); i+=13;
    setbitu(rtcm->buff,i, 4,eph->sva ); i+= 4;
    setbits(rtcm->buff,i,14,idot     ); i+=14;
    setbitu(rtcm->buff,i, 5,eph->iode); i+= 5;
    setbitu(rtcm->buff,i,17,toc      ); i+=17;
    setbits(rtcm->buff,i,11,af2      ); i+=11;
    setbits(rtcm->buff,i,22,af1      ); i+=22;
    setbits(rtcm->buff,i,24,af0      ); i+=24;
    setbitu(rtcm->buff,i, 5,eph->iodc); i+= 5;
    setbits(rtcm->buff,i,18,crs      ); i+=18;
    setbits(rtcm->buff,i,16,deln     ); i+=16;
    setbits(rtcm->buff,i,32,M0       ); i+=32;
    setbits(rtcm->buff,i,18,cuc      ); i+=18;
    setbitu(rtcm->buff,i,32,e        ); i+=32;
    setbits(rtcm->buff,i,18,cus      ); i+=18;
    setbitu(rtcm->buff,i,32,sqrtA    ); i+=32;
    setbitu(rtcm->buff,i,17,toe      ); i+=17;
    setbits(rtcm->buff,i,18,cic      ); i+=18;
    setbits(rtcm->buff,i,32,OMG0     ); i+=32;
    setbits(rtcm->buff,i,18,cis      ); i+=18;
    setbits(rtcm->buff,i,32,i0       ); i+=32;
    setbits(rtcm->buff,i,18,crc      ); i+=18;
    setbits(rtcm->buff,i,32,omg      ); i+=32;
    setbits(rtcm->buff,i,24,OMGd     ); i+=24;
    setbits(rtcm->buff,i,10,tgd1     ); i+=10;
    setbits(rtcm->buff,i,10,tgd2     ); i+=10;
    setbitu(rtcm->buff,i, 1,eph->svh ); i+= 1;
    rtcm->nbit=i;
    return 1;
}
/* encode SSR header ---------------------------------------------------------*/
static int od_rtk_rtcm3e_encode_ssr_head(int type, rtcm_t *rtcm, int sys, int subtype,
                           int nsat, int sync, int iod, double udint, int refd,
                           int provid, int solid)
{
    double tow;
    int i=24,msgno,epoch,week,udi,ns;
    
    trace(4,"encode_ssr_head: type=%d sys=%d subtype=%d nsat=%d sync=%d iod=%d "
          "udint=%.0f\n",type,sys,subtype,nsat,sync,iod,udint);
    
    if (subtype==0) { /* RTCM SSR */
        ns=(sys==SYS_QZS)?4:6;
        switch (sys) {
            case SYS_GPS: msgno=(type==7)?11:1056+type; break;
            case SYS_GLO: msgno=(type==7)? 0:1062+type; break;
            case SYS_GAL: msgno=(type==7)?12:1239+type; break; /* draft */
            case SYS_QZS: msgno=(type==7)?13:1245+type; break; /* draft */
            case SYS_CMP: msgno=(type==7)?14:1257+type; break; /* draft */
            case SYS_SBS: msgno=(type==7)? 0:1251+type; break; /* draft */
            default: return 0;
        }
        if (msgno==0) {
            return 0;
        }
        setbitu(rtcm->buff,i,12,msgno); i+=12; /* message type */
        
        if (sys==SYS_GLO) {
            tow=time2gpst(timeadd(gpst2utc(rtcm->time),10800.0),&week);
            epoch=ROUND(tow)%86400;
            setbitu(rtcm->buff,i,17,epoch); i+=17; /* GLONASS epoch time */
        }
        else {
            tow=time2gpst(rtcm->time,&week);
            epoch=ROUND(tow)%604800;
            setbitu(rtcm->buff,i,20,epoch); i+=20; /* GPS epoch time */
        }
    }
    else { /* IGS SSR */
        ns=6;
        tow=time2gpst(rtcm->time,&week);
        epoch=ROUND(tow)%604800;
        setbitu(rtcm->buff,i,12,4076   ); i+=12; /* message type */
        setbitu(rtcm->buff,i, 3,1      ); i+= 3; /* version */
        setbitu(rtcm->buff,i, 8,subtype); i+= 8; /* subtype */
        setbitu(rtcm->buff,i,20,epoch  ); i+=20; /* SSR epoch time */
    }
    for (udi=0;udi<15;udi++) {
        if (od_rtk_rtcm3e_ssrudint[udi]>=udint) break;
    }
    setbitu(rtcm->buff,i, 4,udi    ); i+= 4; /* update interval */
    setbitu(rtcm->buff,i, 1,sync   ); i+= 1; /* multiple message indicator */
    if (subtype==0&&(type==1||type==4)) {
        setbitu(rtcm->buff,i,1,refd); i+= 1; /* satellite ref datum */
    }
    setbitu(rtcm->buff,i, 4,iod    ); i+= 4; /* IOD SSR */
    setbitu(rtcm->buff,i,16,provid ); i+=16; /* provider ID */
    setbitu(rtcm->buff,i, 4,solid  ); i+= 4; /* solution ID */
    if (subtype>0&&(type==1||type==4)) {
        setbitu(rtcm->buff,i,1,refd); i+= 1; /* global/regional CRS indicator */
    }
    if (type==7) {
        setbitu(rtcm->buff,i,1,0); i+=1; /* dispersive bias consistency ind */
        setbitu(rtcm->buff,i,1,0); i+=1; /* MW consistency indicator */
    }
    setbitu(rtcm->buff,i,ns,nsat); i+=ns; /* no of satellites */
    return i;
}
/* SSR signal and tracking mode IDs ------------------------------------------*/
static  const int od_rtk_rtcm3e_codes_gps[32]={
    CODE_L1C,CODE_L1P,CODE_L1W,CODE_L1S,CODE_L1L,CODE_L2C,CODE_L2D,CODE_L2S,
    CODE_L2L,CODE_L2X,CODE_L2P,CODE_L2W,       0,       0,CODE_L5I,CODE_L5Q
};
static const int od_rtk_rtcm3e_codes_glo[32]={
    CODE_L1C,CODE_L1P,CODE_L2C,CODE_L2P,CODE_L4A,CODE_L4B,CODE_L6A,CODE_L6B,
    CODE_L3I,CODE_L3Q
};
static const int od_rtk_rtcm3e_codes_gal[32]={
    CODE_L1A,CODE_L1B,CODE_L1C,       0,       0,CODE_L5I,CODE_L5Q,       0,
    CODE_L7I,CODE_L7Q,       0,CODE_L8I,CODE_L8Q,       0,CODE_L6A,CODE_L6B,
    CODE_L6C
};
static const int od_rtk_rtcm3e_codes_qzs[32]={
    CODE_L1C,CODE_L1S,CODE_L1L,CODE_L2S,CODE_L2L,       0,CODE_L5I,CODE_L5Q,
           0,CODE_L6S,CODE_L6L,       0,       0,       0,       0,       0,
           0,CODE_L6E
};
static const int od_rtk_rtcm3e_codes_bds[32]={
    CODE_L2I,CODE_L2Q,       0,CODE_L6I,CODE_L6Q,       0,CODE_L7I,CODE_L7Q,
           0,CODE_L1D,CODE_L1P,       0,CODE_L5D,CODE_L5P,       0,CODE_L1A,
           0,       0,CODE_L6A
};
static const int od_rtk_rtcm3e_codes_sbs[32]={
    CODE_L1C,CODE_L5I,CODE_L5Q
};
/* encode SSR 1: orbit corrections -------------------------------------------*/
static int od_rtk_rtcm3e_encode_ssr1(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    double udint=0.0;
    int i,j,iod=0,nsat,prn,iode,iodcrc,refd=0,np,ni,nj,offp,deph[3],ddeph[3];
    
    trace(3,"encode_ssr1: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; ni= 8; nj= 0; offp=  0; break;
        case SYS_GLO: np=5; ni= 8; nj= 0; offp=  0; break;
        case SYS_GAL: np=6; ni=10; nj= 0; offp=  0; break;
        case SYS_QZS: np=4; ni= 8; nj= 0; offp=192; break;
        case SYS_CMP: np=6; ni=10; nj=24; offp=  1; break;
        case SYS_SBS: np=6; ni= 9; nj=24; offp=120; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6; ni=8; nj=0;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[0];
        iod  =rtcm->ssr[j].iod[0];
        refd =rtcm->ssr[j].refd;
    }
    /* encode SSR header */
    i=od_rtk_rtcm3e_encode_ssr_head(1,rtcm,sys,subtype,nsat,sync,iod,udint,refd,0,0);
    
    for (j=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        iode=rtcm->ssr[j].iode;      /* SBAS/BDS: toe/t0 modulo */
        iodcrc=rtcm->ssr[j].iodcrc;  /* SBAS/BDS: IOD CRC */
        
        if (subtype>0) { /* IGS SSR */
            iode&=0xFF;
        }
        deph [0]=ROUND(rtcm->ssr[j].deph [0]/1E-4);
        deph [1]=ROUND(rtcm->ssr[j].deph [1]/4E-4);
        deph [2]=ROUND(rtcm->ssr[j].deph [2]/4E-4);
        ddeph[0]=ROUND(rtcm->ssr[j].ddeph[0]/1E-6);
        ddeph[1]=ROUND(rtcm->ssr[j].ddeph[1]/4E-6);
        ddeph[2]=ROUND(rtcm->ssr[j].ddeph[2]/4E-6);
        
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite ID */
        setbitu(rtcm->buff,i,ni,iode    ); i+=ni; /* IODE */
        setbitu(rtcm->buff,i,nj,iodcrc  ); i+=nj; /* IODCRC */
        setbits(rtcm->buff,i,22,deph [0]); i+=22; /* delta radial */
        setbits(rtcm->buff,i,20,deph [1]); i+=20; /* delta along-track */
        setbits(rtcm->buff,i,20,deph [2]); i+=20; /* delta cross-track */
        setbits(rtcm->buff,i,21,ddeph[0]); i+=21; /* dot delta radial */
        setbits(rtcm->buff,i,19,ddeph[1]); i+=19; /* dot delta along-track */
        setbits(rtcm->buff,i,19,ddeph[2]); i+=19; /* dot delta cross-track */
    }
    rtcm->nbit=i;
    return 1;
}
/* encode SSR 2: clock corrections -------------------------------------------*/
static int od_rtk_rtcm3e_encode_ssr2(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    double udint=0.0;
    int i,j,iod=0,nsat,prn,np,offp,dclk[3];
    
    trace(3,"encode_ssr2: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=192; break;
        case SYS_CMP: np=6; offp=  1; break;
        case SYS_SBS: np=6; offp=120; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[1];
        iod  =rtcm->ssr[j].iod[1];
    }
    /* encode SSR header */
    i=od_rtk_rtcm3e_encode_ssr_head(2,rtcm,sys,subtype,nsat,sync,iod,udint,0,0,0);
    
    for (j=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        dclk[0]=ROUND(rtcm->ssr[j].dclk[0]/1E-4);
        dclk[1]=ROUND(rtcm->ssr[j].dclk[1]/1E-6);
        dclk[2]=ROUND(rtcm->ssr[j].dclk[2]/2E-8);
        
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite ID */
        setbits(rtcm->buff,i,22,dclk[0] ); i+=22; /* delta clock C0 */
        setbits(rtcm->buff,i,21,dclk[1] ); i+=21; /* delta clock C1 */
        setbits(rtcm->buff,i,27,dclk[2] ); i+=27; /* delta clock C2 */
    }
    rtcm->nbit=i;
    return 1;
}
/* encode SSR 3: satellite code biases ---------------------------------------*/
static int od_rtk_rtcm3e_encode_ssr3(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    const int *codes;
    double udint=0.0;
    int i,j,k,iod=0,nsat,prn,nbias,np,offp;
    int code[MAXCODE],bias[MAXCODE];
    
    trace(3,"encode_ssr3: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; codes=od_rtk_rtcm3e_codes_gps; break;
        case SYS_GLO: np=5; offp=  0; codes=od_rtk_rtcm3e_codes_glo; break;
        case SYS_GAL: np=6; offp=  0; codes=od_rtk_rtcm3e_codes_gal; break;
        case SYS_QZS: np=4; offp=192; codes=od_rtk_rtcm3e_codes_qzs; break;
        case SYS_CMP: np=6; offp=  1; codes=od_rtk_rtcm3e_codes_bds; break;
        case SYS_SBS: np=6; offp=120; codes=od_rtk_rtcm3e_codes_sbs; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[4];
        iod  =rtcm->ssr[j].iod[4];
    }
    /* encode SSR header */
    i=od_rtk_rtcm3e_encode_ssr_head(3,rtcm,sys,subtype,nsat,sync,iod,udint,0,0,0);
    
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        for (k=nbias=0;k<32;k++) {
            if (!codes[k]||rtcm->ssr[j].cbias[codes[k]-1]==0.0) continue;
            code[nbias]=k;
            bias[nbias++]=ROUND(rtcm->ssr[j].cbias[codes[k]-1]/0.01);
        }
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite ID */
        setbitu(rtcm->buff,i, 5,nbias);    i+= 5; /* number of code biases */
        
        for (k=0;k<nbias;k++) {
            setbitu(rtcm->buff,i, 5,code[k]); i+= 5; /* signal indicator */
            setbits(rtcm->buff,i,14,bias[k]); i+=14; /* code bias */
        }
    }
    rtcm->nbit=i;
    return 1;
}
/* encode SSR 4: combined orbit and clock corrections ------------------------*/
static int od_rtk_rtcm3e_encode_ssr4(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    double udint=0.0;
    int i,j,iod=0,nsat,prn,iode,iodcrc,refd=0,np,ni,nj,offp;
    int deph[3],ddeph[3],dclk[3];
    
    trace(3,"encode_ssr4: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; ni= 8; nj= 0; offp=  0; break;
        case SYS_GLO: np=5; ni= 8; nj= 0; offp=  0; break;
        case SYS_GAL: np=6; ni=10; nj= 0; offp=  0; break;
        case SYS_QZS: np=4; ni= 8; nj= 0; offp=192; break;
        case SYS_CMP: np=6; ni=10; nj=24; offp=  1; break;
        case SYS_SBS: np=6; ni= 9; nj=24; offp=120; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6; ni=8; nj=0;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[0];
        iod  =rtcm->ssr[j].iod[0];
        refd =rtcm->ssr[j].refd;
    }
    /* encode SSR header */
    i=od_rtk_rtcm3e_encode_ssr_head(4,rtcm,sys,subtype,nsat,sync,iod,udint,refd,0,0);
    
    for (j=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        iode=rtcm->ssr[j].iode;
        iodcrc=rtcm->ssr[j].iodcrc;
        
        if (subtype>0) { /* IGS SSR */
            iode&=0xFF;
        }
        deph [0]=ROUND(rtcm->ssr[j].deph [0]/1E-4);
        deph [1]=ROUND(rtcm->ssr[j].deph [1]/4E-4);
        deph [2]=ROUND(rtcm->ssr[j].deph [2]/4E-4);
        ddeph[0]=ROUND(rtcm->ssr[j].ddeph[0]/1E-6);
        ddeph[1]=ROUND(rtcm->ssr[j].ddeph[1]/4E-6);
        ddeph[2]=ROUND(rtcm->ssr[j].ddeph[2]/4E-6);
        dclk [0]=ROUND(rtcm->ssr[j].dclk [0]/1E-4);
        dclk [1]=ROUND(rtcm->ssr[j].dclk [1]/1E-6);
        dclk [2]=ROUND(rtcm->ssr[j].dclk [2]/2E-8);
        
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite ID */
        setbitu(rtcm->buff,i,ni,iode    ); i+=ni; /* IODE */
        setbitu(rtcm->buff,i,nj,iodcrc  ); i+=nj; /* IODCRC */
        setbits(rtcm->buff,i,22,deph [0]); i+=22; /* delta raidal */
        setbits(rtcm->buff,i,20,deph [1]); i+=20; /* delta along-track */
        setbits(rtcm->buff,i,20,deph [2]); i+=20; /* delta cross-track */
        setbits(rtcm->buff,i,21,ddeph[0]); i+=21; /* dot delta radial */
        setbits(rtcm->buff,i,19,ddeph[1]); i+=19; /* dot delta along-track */
        setbits(rtcm->buff,i,19,ddeph[2]); i+=19; /* dot delta cross-track */
        setbits(rtcm->buff,i,22,dclk [0]); i+=22; /* delta clock C0 */
        setbits(rtcm->buff,i,21,dclk [1]); i+=21; /* delta clock C1 */
        setbits(rtcm->buff,i,27,dclk [2]); i+=27; /* delta clock C2 */
    }
    rtcm->nbit=i;
    return 1;
}
/* encode SSR 5: URA ---------------------------------------------------------*/
static int od_rtk_rtcm3e_encode_ssr5(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    double udint=0.0;
    int i,j,nsat,iod=0,prn,ura,np,offp;
    
    trace(3,"encode_ssr5: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=192; break;
        case SYS_CMP: np=6; offp=  1; break;
        case SYS_SBS: np=6; offp=120; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[3];
        iod  =rtcm->ssr[j].iod[3];
    }
    /* encode ssr header */
    i=od_rtk_rtcm3e_encode_ssr_head(5,rtcm,sys,subtype,nsat,sync,iod,udint,0,0,0);
    
    for (j=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        ura=rtcm->ssr[j].ura;
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite id */
        setbitu(rtcm->buff,i, 6,ura     ); i+= 6; /* ssr ura */
    }
    rtcm->nbit=i;
    return 1;
}
/* encode SSR 6: high rate clock correction ----------------------------------*/
static int od_rtk_rtcm3e_encode_ssr6(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    double udint=0.0;
    int i,j,nsat,iod=0,prn,hrclk,np,offp;
    
    trace(3,"encode_ssr6: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; break;
        case SYS_GLO: np=5; offp=  0; break;
        case SYS_GAL: np=6; offp=  0; break;
        case SYS_QZS: np=4; offp=192; break;
        case SYS_CMP: np=6; offp=  1; break;
        case SYS_SBS: np=6; offp=120; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[2];
        iod  =rtcm->ssr[j].iod[2];
    }
    /* encode SSR header */
    i=od_rtk_rtcm3e_encode_ssr_head(6,rtcm,sys,subtype,nsat,sync,iod,udint,0,0,0);
    
    for (j=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        hrclk=ROUND(rtcm->ssr[j].hrclk/1E-4);
        
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite ID */
        setbits(rtcm->buff,i,22,hrclk   ); i+=22; /* high rate clock corr */
    }
    rtcm->nbit=i;
    return 1;
}
/* encode SSR 7: satellite phase biases --------------------------------------*/
static int od_rtk_rtcm3e_encode_ssr7(rtcm_t *rtcm, int sys, int subtype, int sync)
{
    const int *codes;
    double udint=0.0;
    int i,j,k,iod=0,nsat,prn,nbias,np,offp;
    int code[MAXCODE],pbias[MAXCODE],stdpb[MAXCODE],yaw_ang,yaw_rate;
    
    trace(3,"encode_ssr7: sys=%d subtype=%d sync=%d\n",sys,subtype,sync);
    
    switch (sys) {
        case SYS_GPS: np=6; offp=  0; codes=od_rtk_rtcm3e_codes_gps; break;
        case SYS_GLO: np=5; offp=  0; codes=od_rtk_rtcm3e_codes_glo; break;
        case SYS_GAL: np=6; offp=  0; codes=od_rtk_rtcm3e_codes_gal; break;
        case SYS_QZS: np=4; offp=192; codes=od_rtk_rtcm3e_codes_qzs; break;
        case SYS_CMP: np=6; offp=  1; codes=od_rtk_rtcm3e_codes_bds; break;
        case SYS_SBS: np=6; offp=120; codes=od_rtk_rtcm3e_codes_sbs; break;
        default: return 0;
    }
    if (subtype>0) { /* IGS SSR */
        np=6;
        if      (sys==SYS_CMP) offp=0;
        else if (sys==SYS_SBS) offp=119;
    }
    /* number of satellites */
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        nsat++;
        udint=rtcm->ssr[j].udi[5];
        iod  =rtcm->ssr[j].iod[5];
    }
    /* encode SSR header */
    i=od_rtk_rtcm3e_encode_ssr_head(7,rtcm,sys,subtype,nsat,sync,iod,udint,0,0,0);
    
    for (j=nsat=0;j<MAXSAT;j++) {
        if (satsys(j+1,&prn)!=sys||!rtcm->ssr[j].update) continue;
        
        for (k=nbias=0;k<32;k++) {
            if (!codes[k]||rtcm->ssr[j].pbias[codes[k]-1]==0.0) continue;
            code[nbias]=k;
            pbias[nbias  ]=ROUND(rtcm->ssr[j].pbias[codes[k]-1]/0.0001);
            stdpb[nbias++]=ROUND(rtcm->ssr[j].stdpb[codes[k]-1]/0.0001);
        }
        yaw_ang =ROUND(rtcm->ssr[j].yaw_ang /180.0* 256.0);
        yaw_rate=ROUND(rtcm->ssr[j].yaw_rate/180.0*8192.0);
        setbitu(rtcm->buff,i,np,prn-offp); i+=np; /* satellite ID */
        setbitu(rtcm->buff,i, 5,nbias);    i+= 5; /* number of code biases */
        setbitu(rtcm->buff,i, 9,yaw_ang);  i+= 9; /* yaw angle */
        setbits(rtcm->buff,i, 8,yaw_rate); i+= 8; /* yaw rate */
        
        for (k=0;k<nbias;k++) {
            setbitu(rtcm->buff,i, 5,code[k] ); i+= 5; /* signal indicator */
            setbitu(rtcm->buff,i, 1,0       ); i+= 1; /* integer-indicator */
            setbitu(rtcm->buff,i, 2,0       ); i+= 2; /* WL integer-indicator */
            setbitu(rtcm->buff,i, 4,0       ); i+= 4; /* discont counter */
            setbits(rtcm->buff,i,20,pbias[k]); i+=20; /* phase bias */
            if (subtype==0) {
                setbits(rtcm->buff,i,17,stdpb[k]); i+=17; /* std-dev ph-bias */
            }
        }
    }
    rtcm->nbit=i;
    return 1;
}
/* satellite no to MSM satellite ID ------------------------------------------*/
static int od_rtk_rtcm3e_to_satid(int sys, int sat)
{
    int prn;
    
    if (satsys(sat,&prn)!=sys) return 0;
    
    if      (sys==SYS_QZS) prn-=MINPRNQZS-1;
    else if (sys==SYS_SBS) prn-=MINPRNSBS-1;
    
    return prn;
}
/* observation code to MSM signal ID -----------------------------------------*/
static int od_rtk_rtcm3e_to_sigid(int sys, uint8_t code)
{
    const char **msm_sig;
    char *sig;
    int i;
    
    /* signal conversion for undefined signal by rtcm */
    if (sys==SYS_GPS) {
        if      (code==CODE_L1Y) code=CODE_L1P;
        else if (code==CODE_L1M) code=CODE_L1P;
        else if (code==CODE_L1N) code=CODE_L1P;
        else if (code==CODE_L2D) code=CODE_L2P;
        else if (code==CODE_L2Y) code=CODE_L2P;
        else if (code==CODE_L2M) code=CODE_L2P;
        else if (code==CODE_L2N) code=CODE_L2P;
    }
    if (!*(sig=code2obs(code))) return 0;
    
    switch (sys) {
        case SYS_GPS: msm_sig=msm_sig_gps; break;
        case SYS_GLO: msm_sig=msm_sig_glo; break;
        case SYS_GAL: msm_sig=msm_sig_gal; break;
        case SYS_QZS: msm_sig=msm_sig_qzs; break;
        case SYS_SBS: msm_sig=msm_sig_sbs; break;
        case SYS_CMP: msm_sig=msm_sig_cmp; break;
        case SYS_IRN: msm_sig=msm_sig_irn; break;
        default: return 0;
    }
    for (i=0;i<32;i++) {
        if (!strcmp(sig,msm_sig[i])) return i+1;
    }
    return 0;
}
/* generate MSM satellite, signal and cell index -----------------------------*/
static void od_rtk_rtcm3e_gen_msm_index(rtcm_t *rtcm, int sys, int *nsat, int *nsig,
                          int *ncell, uint8_t *sat_ind, uint8_t *sig_ind,
                          uint8_t *cell_ind)
{
    int i,j,sat,sig,cell;
    
    *nsat=*nsig=*ncell=0;
    
    /* generate satellite and signal index */
    for (i=0;i<rtcm->obs.n;i++) {
        if (!(sat=od_rtk_rtcm3e_to_satid(sys,rtcm->obs.data[i].sat))) continue;
        
        for (j=0;j<NFREQ+NEXOBS;j++) {
            if (!(sig=od_rtk_rtcm3e_to_sigid(sys,rtcm->obs.data[i].code[j]))) continue;
            
            sat_ind[sat-1]=sig_ind[sig-1]=1;
        }
    }
    for (i=0;i<64;i++) {
        if (sat_ind[i]) sat_ind[i]=++(*nsat);
    }
    for (i=0;i<32;i++) {
        if (sig_ind[i]) sig_ind[i]=++(*nsig);
    }
    /* generate cell index */
    for (i=0;i<rtcm->obs.n;i++) {
        if (!(sat=od_rtk_rtcm3e_to_satid(sys,rtcm->obs.data[i].sat))) continue;
        
        for (j=0;j<NFREQ+NEXOBS;j++) {
            if (!(sig=od_rtk_rtcm3e_to_sigid(sys,rtcm->obs.data[i].code[j]))) continue;
            
            cell=sig_ind[sig-1]-1+(sat_ind[sat-1]-1)*(*nsig);
            cell_ind[cell]=1;
        }
    }
    for (i=0;i<*nsat*(*nsig);i++) {
        if (cell_ind[i]&&*ncell<64) cell_ind[i]=++(*ncell);
    }
}
/* generate MSM satellite data fields ----------------------------------------*/
static void od_rtk_rtcm3e_gen_msm_sat(rtcm_t *rtcm, int sys, int nsat, const uint8_t *sat_ind,
                        double *rrng, double *rrate, uint8_t *info)
{
    obsd_t *data;
    double freq;
    int i,j,k,sat,sig,fcn;
    
    for (i=0;i<64;i++) rrng[i]=rrate[i]=0.0;
    
    for (i=0;i<rtcm->obs.n;i++) {
        data=rtcm->obs.data+i;
        fcn=od_rtk_rtcm3e_fcn_glo(data->sat,rtcm); /* fcn+7 */
        
        if (!(sat=od_rtk_rtcm3e_to_satid(sys,data->sat))) continue;
        
        for (j=0;j<NFREQ+NEXOBS;j++) {
            if (!(sig=od_rtk_rtcm3e_to_sigid(sys,data->code[j]))) continue;
            k=sat_ind[sat-1]-1;
            freq=code2freq(sys,data->code[j],fcn-7);
            
            /* rough range (ms) and rough phase-range-rate (m/s) */
            if (rrng[k]==0.0&&data->P[j]!=0.0) {
                rrng[k]=ROUND( data->P[j]/RANGE_MS/P2_10)*RANGE_MS*P2_10;
            }
            if (rrate[k]==0.0&&data->D[j]!=0.0&&freq>0.0) {
                rrate[k]=ROUND(-data->D[j]*CLIGHT/freq)*1.0;
            }
            /* extended satellite info */
            if (info) info[k]=sys!=SYS_GLO?0:(fcn<0?15:fcn);
        }
    }
}
/* generate MSM signal data fields -------------------------------------------*/
static void od_rtk_rtcm3e_gen_msm_sig(rtcm_t *rtcm, int sys, int nsat, int nsig, int ncell,
                        const uint8_t *sat_ind, const uint8_t *sig_ind,
                        const uint8_t *cell_ind, const double *rrng,
                        const double *rrate, double *psrng, double *phrng,
                        double *rate, double *lock, uint8_t *half, float *cnr)
{
    obsd_t *data;
    double freq,lambda,psrng_s,phrng_s,rate_s,lt;
    int i,j,k,sat,sig,fcn,cell,LLI;
    
    for (i=0;i<ncell;i++) {
        if (psrng) psrng[i]=0.0;
        if (phrng) phrng[i]=0.0;
        if (rate ) rate [i]=0.0;
    }
    for (i=0;i<rtcm->obs.n;i++) {
        data=rtcm->obs.data+i;
        fcn=od_rtk_rtcm3e_fcn_glo(data->sat,rtcm); /* fcn+7 */
        
        if (!(sat=od_rtk_rtcm3e_to_satid(sys,data->sat))) continue;
        
        for (j=0;j<NFREQ+NEXOBS;j++) {
            if (!(sig=od_rtk_rtcm3e_to_sigid(sys,data->code[j]))) continue;
            
            k=sat_ind[sat-1]-1;
            if ((cell=cell_ind[sig_ind[sig-1]-1+k*nsig])>=64) continue;
            
            freq=code2freq(sys,data->code[j],fcn-7);
            lambda=freq==0.0?0.0:CLIGHT/freq;
            psrng_s=data->P[j]==0.0?0.0:data->P[j]-rrng[k];
            phrng_s=data->L[j]==0.0||lambda<=0.0?0.0: data->L[j]*lambda-rrng [k];
            rate_s =data->D[j]==0.0||lambda<=0.0?0.0:-data->D[j]*lambda-rrate[k];
            
            /* subtract phase - psudorange integer cycle offset */
            LLI=data->LLI[j];
            if ((LLI&1)||fabs(phrng_s-rtcm->cp[data->sat-1][j])>1171.0) {
                rtcm->cp[data->sat-1][j]=ROUND(phrng_s/lambda)*lambda;
                LLI|=1;
            }
            phrng_s-=rtcm->cp[data->sat-1][j];
            
            lt=od_rtk_rtcm3e_locktime_d(data->time,rtcm->lltime[data->sat-1]+j,LLI);
            
            if (psrng&&psrng_s!=0.0) psrng[cell-1]=psrng_s;
            if (phrng&&phrng_s!=0.0) phrng[cell-1]=phrng_s;
            if (rate &&rate_s !=0.0) rate [cell-1]=rate_s;
            if (lock) lock[cell-1]=lt;
            if (half) half[cell-1]=(data->LLI[j]&2)?1:0;
            if (cnr ) cnr [cell-1]=(float)(data->SNR[j]*SNR_UNIT);
        }
    }
}
/* encode MSM header ---------------------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_head(int type, rtcm_t *rtcm, int sys, int sync, int *nsat,
                           int *ncell, double *rrng, double *rrate,
                           uint8_t *info, double *psrng, double *phrng,
                           double *rate, double *lock, uint8_t *half,
                           float *cnr)
{
    double tow;
    uint8_t sat_ind[64]={0},sig_ind[32]={0},cell_ind[32*64]={0};
    uint32_t dow,epoch;
    int i=24,j,nsig=0;
    
    switch (sys) {
        case SYS_GPS: type+=1070; break;
        case SYS_GLO: type+=1080; break;
        case SYS_GAL: type+=1090; break;
        case SYS_QZS: type+=1110; break;
        case SYS_SBS: type+=1100; break;
        case SYS_CMP: type+=1120; break;
        case SYS_IRN: type+=1130; break;
        default: return 0;
    }
    /* generate msm satellite, signal and cell index */
    od_rtk_rtcm3e_gen_msm_index(rtcm,sys,nsat,&nsig,ncell,sat_ind,sig_ind,cell_ind);
    
    if (sys==SYS_GLO) {
        /* GLONASS time (dow + tod-ms) */
        tow=time2gpst(timeadd(gpst2utc(rtcm->time),10800.0),NULL);
        dow=(uint32_t)(tow/86400.0);
        epoch=(dow<<27)+ROUND_U(fmod(tow,86400.0)*1E3);
    }
    else if (sys==SYS_CMP) {
        /* BDS time (tow-ms) */
        epoch=ROUND_U(time2gpst(gpst2bdt(rtcm->time),NULL)*1E3);
    }
    else {
        /* GPS, QZSS, Galileo and IRNSS time (tow-ms) */
        epoch=ROUND_U(time2gpst(rtcm->time,NULL)*1E3);
    }
    /* encode msm header (ref [15] table 3.5-78) */
    setbitu(rtcm->buff,i,12,type       ); i+=12; /* message number */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* reference station id */
    setbitu(rtcm->buff,i,30,epoch      ); i+=30; /* epoch time */
    setbitu(rtcm->buff,i, 1,sync       ); i+= 1; /* multiple message bit */
    setbitu(rtcm->buff,i, 3,rtcm->seqno); i+= 3; /* issue of data station */
    setbitu(rtcm->buff,i, 7,0          ); i+= 7; /* reserved */
    setbitu(rtcm->buff,i, 2,0          ); i+= 2; /* clock streering indicator */
    setbitu(rtcm->buff,i, 2,0          ); i+= 2; /* external clock indicator */
    setbitu(rtcm->buff,i, 1,0          ); i+= 1; /* smoothing indicator */
    setbitu(rtcm->buff,i, 3,0          ); i+= 3; /* smoothing interval */
    
    /* satellite mask */
    for (j=0;j<64;j++) {
        setbitu(rtcm->buff,i,1,sat_ind[j]?1:0); i+=1;
    }
    /* signal mask */
    for (j=0;j<32;j++) {
        setbitu(rtcm->buff,i,1,sig_ind[j]?1:0); i+=1;
    }
    /* cell mask */
    for (j=0;j<*nsat*nsig&&j<64;j++) {
        setbitu(rtcm->buff,i,1,cell_ind[j]?1:0); i+=1;
    }
    /* generate msm satellite data fields */
    od_rtk_rtcm3e_gen_msm_sat(rtcm,sys,*nsat,sat_ind,rrng,rrate,info);
    
    /* generate msm signal data fields */
    od_rtk_rtcm3e_gen_msm_sig(rtcm,sys,*nsat,nsig,*ncell,sat_ind,sig_ind,cell_ind,rrng,rrate,
                psrng,phrng,rate,lock,half,cnr);
    
    return i;
}
/* encode rough range integer ms ---------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_int_rrng(rtcm_t *rtcm, int i, const double *rrng,
                               int nsat)
{
    uint32_t int_ms;
    int j;
    
    for (j=0;j<nsat;j++) {
        if (rrng[j]==0.0) {
            int_ms=255;
        }
        else if (rrng[j]<0.0||rrng[j]>RANGE_MS*255.0) {
            trace(2,"msm rough range overflow %s rrng=%.3f\n",
                 time_str(rtcm->time,0),rrng[j]);
            int_ms=255;
        }
        else {
            int_ms=ROUND_U(rrng[j]/RANGE_MS/P2_10)>>10;
        }
        setbitu(rtcm->buff,i,8,int_ms); i+=8;
    }
    return i;
}
/* encode rough range modulo 1 ms --------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm_t *rtcm, int i, const double *rrng,
                               int nsat)
{
    uint32_t mod_ms;
    int j;
    
    for (j=0;j<nsat;j++) {
        if (rrng[j]<=0.0||rrng[j]>RANGE_MS*255.0) {
            mod_ms=0;
        }
        else {
            mod_ms=ROUND_U(rrng[j]/RANGE_MS/P2_10)&0x3FFu;
        }
        setbitu(rtcm->buff,i,10,mod_ms); i+=10;
    }
    return i;
}
/* encode extended satellite info --------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_info(rtcm_t *rtcm, int i, const uint8_t *info, int nsat)
{
    int j;
    
    for (j=0;j<nsat;j++) {
        setbitu(rtcm->buff,i,4,info[j]); i+=4;
    }
    return i;
}
/* encode rough phase-range-rate ---------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_rrate(rtcm_t *rtcm, int i, const double *rrate, int nsat)
{
    int j,rrate_val;
    
    for (j=0;j<nsat;j++) {
        if (fabs(rrate[j])>8191.0) {
            trace(2,"msm rough phase-range-rate overflow %s rrate=%.4f\n",
                 time_str(rtcm->time,0),rrate[j]);
            rrate_val=-8192;
        }
        else {
            rrate_val=ROUND(rrate[j]/1.0);
        }
        setbits(rtcm->buff,i,14,rrate_val); i+=14;
    }
    return i;
}
/* encode fine pseudorange ---------------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_psrng(rtcm_t *rtcm, int i, const double *psrng, int ncell)
{
    int j,psrng_val;
    
    for (j=0;j<ncell;j++) {
        if (psrng[j]==0.0) {
            psrng_val=-16384;
        }
        else if (fabs(psrng[j])>292.7) {
            trace(2,"msm fine pseudorange overflow %s psrng=%.3f\n",
                 time_str(rtcm->time,0),psrng[j]);
            psrng_val=-16384;
        }
        else {
            psrng_val=ROUND(psrng[j]/RANGE_MS/P2_24);
        }
        setbits(rtcm->buff,i,15,psrng_val); i+=15;
    }
    return i;
}
/* encode fine pseudorange with extended resolution --------------------------*/
static int od_rtk_rtcm3e_encode_msm_psrng_ex(rtcm_t *rtcm, int i, const double *psrng,
                               int ncell)
{
    int j,psrng_val;
    
    for (j=0;j<ncell;j++) {
        if (psrng[j]==0.0) {
            psrng_val=-524288;
        }
        else if (fabs(psrng[j])>292.7) {
            trace(2,"msm fine pseudorange ext overflow %s psrng=%.3f\n",
                 time_str(rtcm->time,0),psrng[j]);
            psrng_val=-524288;
        }
        else {
            psrng_val=ROUND(psrng[j]/RANGE_MS/P2_29);
        }
        setbits(rtcm->buff,i,20,psrng_val); i+=20;
    }
    return i;
}
/* encode fine phase-range ---------------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_phrng(rtcm_t *rtcm, int i, const double *phrng, int ncell)
{
    int j,phrng_val;
    
    for (j=0;j<ncell;j++) {
        if (phrng[j]==0.0) {
            phrng_val=-2097152;
        }
        else if (fabs(phrng[j])>1171.0) {
            trace(2,"msm fine phase-range overflow %s phrng=%.3f\n",
                 time_str(rtcm->time,0),phrng[j]);
            phrng_val=-2097152;
        }
        else {
            phrng_val=ROUND(phrng[j]/RANGE_MS/P2_29);
        }
        setbits(rtcm->buff,i,22,phrng_val); i+=22;
    }
    return i;
}
/* encode fine phase-range with extended resolution --------------------------*/
static int od_rtk_rtcm3e_encode_msm_phrng_ex(rtcm_t *rtcm, int i, const double *phrng,
                               int ncell)
{
    int j,phrng_val;
    
    for (j=0;j<ncell;j++) {
        if (phrng[j]==0.0) {
            phrng_val=-8388608;
        }
        else if (fabs(phrng[j])>1171.0) {
            trace(2,"msm fine phase-range ext overflow %s phrng=%.3f\n",
                 time_str(rtcm->time,0),phrng[j]);
            phrng_val=-8388608;
        }
        else {
            phrng_val=ROUND(phrng[j]/RANGE_MS/P2_31);
        }
        setbits(rtcm->buff,i,24,phrng_val); i+=24;
    }
    return i;
}
/* encode lock-time indicator ------------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_lock(rtcm_t *rtcm, int i, const double *lock, int ncell)
{
    int j,lock_val;
    
    for (j=0;j<ncell;j++) {
        lock_val=od_rtk_rtcm3e_to_msm_lock(lock[j]);
        setbitu(rtcm->buff,i,4,lock_val); i+=4;
    }
    return i;
}
/* encode lock-time indicator with extended range and resolution -------------*/
static int od_rtk_rtcm3e_encode_msm_lock_ex(rtcm_t *rtcm, int i, const double *lock,
                              int ncell)
{
    int j,lock_val;
    
    for (j=0;j<ncell;j++) {
        lock_val=od_rtk_rtcm3e_to_msm_lock_ex(lock[j]);
        setbitu(rtcm->buff,i,10,lock_val); i+=10;
    }
    return i;
}
/* encode half-cycle-ambiguity indicator -------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_half_amb(rtcm_t *rtcm, int i, const uint8_t *half,
                               int ncell)
{
    int j;
    
    for (j=0;j<ncell;j++) {
        setbitu(rtcm->buff,i,1,half[j]); i+=1;
    }
    return i;
}
/* encode signal CNR ---------------------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_cnr(rtcm_t *rtcm, int i, const float *cnr, int ncell)
{
    int j,cnr_val;
    
    for (j=0;j<ncell;j++) {
        cnr_val=ROUND(cnr[j]/1.0);
        setbitu(rtcm->buff,i,6,cnr_val); i+=6;
    }
    return i;
}
/* encode signal CNR with extended resolution --------------------------------*/
static int od_rtk_rtcm3e_encode_msm_cnr_ex(rtcm_t *rtcm, int i, const float *cnr, int ncell)
{
    int j,cnr_val;
    
    for (j=0;j<ncell;j++) {
        cnr_val=ROUND(cnr[j]/0.0625);
        setbitu(rtcm->buff,i,10,cnr_val); i+=10;
    }
    return i;
}
/* encode fine phase-range-rate ----------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm_rate(rtcm_t *rtcm, int i, const double *rate, int ncell)
{
    int j,rate_val;
    
    for (j=0;j<ncell;j++) {
        if (rate[j]==0.0) {
            rate_val=-16384;
        }
        else if (fabs(rate[j])>1.6384) {
            trace(2,"msm fine phase-range-rate overflow %s rate=%.3f\n",
                 time_str(rtcm->time,0),rate[j]);
            rate_val=-16384;
        }
        else {
            rate_val=ROUND(rate[j]/0.0001);
        }
        setbitu(rtcm->buff,i,15,rate_val); i+=15;
    }
    return i;
}
/* encode MSM 1: compact pseudorange -----------------------------------------*/
static int od_rtk_rtcm3e_encode_msm1(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],psrng[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm1: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(1,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,NULL,psrng,
                            NULL,NULL,NULL,NULL,NULL))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_psrng   (rtcm,i,psrng,ncell); /* fine pseudorange */
    
    rtcm->nbit=i;
    return 1;
}
/* encode MSM 2: compact phaserange ------------------------------------------*/
static int od_rtk_rtcm3e_encode_msm2(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],phrng[64],lock[64];
    uint8_t half[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm2: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(2,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,NULL,NULL,
                            phrng,NULL,lock,half,NULL))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_phrng   (rtcm,i,phrng,ncell); /* fine phase-range */
    i=od_rtk_rtcm3e_encode_msm_lock    (rtcm,i,lock ,ncell); /* lock-time indicator */
    i=od_rtk_rtcm3e_encode_msm_half_amb(rtcm,i,half ,ncell); /* half-cycle-amb indicator */
    
    rtcm->nbit=i;
    return 1;
}
/* encode MSM 3: compact pseudorange and phaserange --------------------------*/
static int od_rtk_rtcm3e_encode_msm3(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],psrng[64],phrng[64],lock[64];
    uint8_t half[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm3: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(3,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,NULL,psrng,
                            phrng,NULL,lock,half,NULL))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_psrng   (rtcm,i,psrng,ncell); /* fine pseudorange */
    i=od_rtk_rtcm3e_encode_msm_phrng   (rtcm,i,phrng,ncell); /* fine phase-range */
    i=od_rtk_rtcm3e_encode_msm_lock    (rtcm,i,lock ,ncell); /* lock-time indicator */
    i=od_rtk_rtcm3e_encode_msm_half_amb(rtcm,i,half ,ncell); /* half-cycle-amb indicator */
    
    rtcm->nbit=i;
    return 1;
}
/* encode MSM 4: full pseudorange and phaserange plus CNR --------------------*/
static int od_rtk_rtcm3e_encode_msm4(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],psrng[64],phrng[64],lock[64];
    float cnr[64];
    uint8_t half[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm4: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(4,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,NULL,psrng,
                            phrng,NULL,lock,half,cnr))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_int_rrng(rtcm,i,rrng ,nsat ); /* rough range integer ms */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_psrng   (rtcm,i,psrng,ncell); /* fine pseudorange */
    i=od_rtk_rtcm3e_encode_msm_phrng   (rtcm,i,phrng,ncell); /* fine phase-range */
    i=od_rtk_rtcm3e_encode_msm_lock    (rtcm,i,lock ,ncell); /* lock-time indicator */
    i=od_rtk_rtcm3e_encode_msm_half_amb(rtcm,i,half ,ncell); /* half-cycle-amb indicator */
    i=od_rtk_rtcm3e_encode_msm_cnr     (rtcm,i,cnr  ,ncell); /* signal cnr */
    rtcm->nbit=i;
    return 1;
}
/* encode MSM 5: full pseudorange, phaserange, phaserangerate and CNR --------*/
static int od_rtk_rtcm3e_encode_msm5(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],psrng[64],phrng[64],rate[64],lock[64];
    float cnr[64];
    uint8_t info[64],half[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm5: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(5,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,info,psrng,
                            phrng,rate,lock,half,cnr))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_int_rrng(rtcm,i,rrng ,nsat ); /* rough range integer ms */
    i=od_rtk_rtcm3e_encode_msm_info    (rtcm,i,info ,nsat ); /* extended satellite info */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    i=od_rtk_rtcm3e_encode_msm_rrate   (rtcm,i,rrate,nsat ); /* rough phase-range-rate */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_psrng   (rtcm,i,psrng,ncell); /* fine pseudorange */
    i=od_rtk_rtcm3e_encode_msm_phrng   (rtcm,i,phrng,ncell); /* fine phase-range */
    i=od_rtk_rtcm3e_encode_msm_lock    (rtcm,i,lock ,ncell); /* lock-time indicator */
    i=od_rtk_rtcm3e_encode_msm_half_amb(rtcm,i,half ,ncell); /* half-cycle-amb indicator */
    i=od_rtk_rtcm3e_encode_msm_cnr     (rtcm,i,cnr  ,ncell); /* signal cnr */
    i=od_rtk_rtcm3e_encode_msm_rate    (rtcm,i,rate ,ncell); /* fine phase-range-rate */
    rtcm->nbit=i;
    return 1;
}
/* encode MSM 6: full pseudorange and phaserange plus CNR (high-res) ---------*/
static int od_rtk_rtcm3e_encode_msm6(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],psrng[64],phrng[64],lock[64];
    float cnr[64];
    uint8_t half[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm6: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(6,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,NULL,psrng,
                            phrng,NULL,lock,half,cnr))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_int_rrng(rtcm,i,rrng ,nsat ); /* rough range integer ms */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_psrng_ex(rtcm,i,psrng,ncell); /* fine pseudorange ext */
    i=od_rtk_rtcm3e_encode_msm_phrng_ex(rtcm,i,phrng,ncell); /* fine phase-range ext */
    i=od_rtk_rtcm3e_encode_msm_lock_ex (rtcm,i,lock ,ncell); /* lock-time indicator ext */
    i=od_rtk_rtcm3e_encode_msm_half_amb(rtcm,i,half ,ncell); /* half-cycle-amb indicator */
    i=od_rtk_rtcm3e_encode_msm_cnr_ex  (rtcm,i,cnr  ,ncell); /* signal cnr ext */
    rtcm->nbit=i;
    return 1;
}
/* encode MSM 7: full pseudorange, phaserange, phaserangerate and CNR (h-res) */
static int od_rtk_rtcm3e_encode_msm7(rtcm_t *rtcm, int sys, int sync)
{
    double rrng[64],rrate[64],psrng[64],phrng[64],rate[64],lock[64];
    float cnr[64];
    uint8_t info[64],half[64];
    int i,nsat,ncell;
    
    trace(3,"encode_msm7: sys=%d sync=%d\n",sys,sync);
    
    /* encode msm header */
    if (!(i=od_rtk_rtcm3e_encode_msm_head(7,rtcm,sys,sync,&nsat,&ncell,rrng,rrate,info,psrng,
                            phrng,rate,lock,half,cnr))) {
        return 0;
    }
    /* encode msm satellite data */
    i=od_rtk_rtcm3e_encode_msm_int_rrng(rtcm,i,rrng ,nsat ); /* rough range integer ms */
    i=od_rtk_rtcm3e_encode_msm_info    (rtcm,i,info ,nsat ); /* extended satellite info */
    i=od_rtk_rtcm3e_encode_msm_mod_rrng(rtcm,i,rrng ,nsat ); /* rough range modulo 1 ms */
    i=od_rtk_rtcm3e_encode_msm_rrate   (rtcm,i,rrate,nsat ); /* rough phase-range-rate */
    
    /* encode msm signal data */
    i=od_rtk_rtcm3e_encode_msm_psrng_ex(rtcm,i,psrng,ncell); /* fine pseudorange ext */
    i=od_rtk_rtcm3e_encode_msm_phrng_ex(rtcm,i,phrng,ncell); /* fine phase-range ext */
    i=od_rtk_rtcm3e_encode_msm_lock_ex (rtcm,i,lock ,ncell); /* lock-time indicator ext */
    i=od_rtk_rtcm3e_encode_msm_half_amb(rtcm,i,half ,ncell); /* half-cycle-amb indicator */
    i=od_rtk_rtcm3e_encode_msm_cnr_ex  (rtcm,i,cnr  ,ncell); /* signal cnr ext */
    i=od_rtk_rtcm3e_encode_msm_rate    (rtcm,i,rate ,ncell); /* fine phase-range-rate */
    rtcm->nbit=i;
    return 1;
}
/* encode type 1230: GLONASS L1 and L2 code-phase biases ---------------------*/
static int od_rtk_rtcm3e_encode_type1230(rtcm_t *rtcm, int sync)
{
    int i=24,j,align,mask=15,bias[4];
    
    trace(3,"encode_type1230: sync=%d\n",sync);
    
    align=rtcm->sta.glo_cp_align;
    
    for (j=0;j<4;j++) {
        bias[j]=ROUND(rtcm->sta.glo_cp_bias[j]/0.02);
        if (bias[j]<=-32768||bias[j]>32767) {
            bias[j]=-32768; /* invalid value */
        }
    }
    setbitu(rtcm->buff,i,12,1230       ); i+=12; /* message no */
    setbitu(rtcm->buff,i,12,rtcm->staid); i+=12; /* station ID */
    setbitu(rtcm->buff,i, 1,align      ); i+= 1; /* GLO code-phase bias ind */
    setbitu(rtcm->buff,i, 3,0          ); i+= 3; /* reserved */
    setbitu(rtcm->buff,i, 4,mask       ); i+= 4; /* GLO FDMA signals mask */
    setbits(rtcm->buff,i,16,bias[0]    ); i+=16; /* GLO C1 code-phase bias */
    setbits(rtcm->buff,i,16,bias[1]    ); i+=16; /* GLO P1 code-phase bias */
    setbits(rtcm->buff,i,16,bias[2]    ); i+=16; /* GLO C2 code-phase bias */
    setbits(rtcm->buff,i,16,bias[3]    ); i+=16; /* GLO P2 code-phase bias */
    rtcm->nbit=i;
    return 1;
}
/* encode type 4073: proprietary message Mitsubishi Electric -----------------*/
static int od_rtk_rtcm3e_encode_type4073(rtcm_t *rtcm, int subtype, int sync)
{
    trace(2,"rtcm3 4073: unsupported message subtype=%d\n",subtype);
    return 0;
}
/* encode type 4076: proprietary message IGS ---------------------------------*/
static int od_rtk_rtcm3e_encode_type4076(rtcm_t *rtcm, int subtype, int sync)
{
    switch (subtype) {
        case  21: return od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_GPS,subtype,sync);
        case  22: return od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_GPS,subtype,sync);
        case  23: return od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_GPS,subtype,sync);
        case  24: return od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_GPS,subtype,sync);
        case  25: return od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_GPS,subtype,sync);
        case  26: return od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_GPS,subtype,sync);
        case  27: return od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_GPS,subtype,sync);
        case  41: return od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_GLO,subtype,sync);
        case  42: return od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_GLO,subtype,sync);
        case  43: return od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_GLO,subtype,sync);
        case  44: return od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_GLO,subtype,sync);
        case  45: return od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_GLO,subtype,sync);
        case  46: return od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_GLO,subtype,sync);
        case  47: return od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_GLO,subtype,sync);
        case  61: return od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_GAL,subtype,sync);
        case  62: return od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_GAL,subtype,sync);
        case  63: return od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_GAL,subtype,sync);
        case  64: return od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_GAL,subtype,sync);
        case  65: return od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_GAL,subtype,sync);
        case  66: return od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_GAL,subtype,sync);
        case  67: return od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_GAL,subtype,sync);
        case  81: return od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_QZS,subtype,sync);
        case  82: return od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_QZS,subtype,sync);
        case  83: return od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_QZS,subtype,sync);
        case  84: return od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_QZS,subtype,sync);
        case  85: return od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_QZS,subtype,sync);
        case  86: return od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_QZS,subtype,sync);
        case  87: return od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_QZS,subtype,sync);
        case 101: return od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_CMP,subtype,sync);
        case 102: return od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_CMP,subtype,sync);
        case 103: return od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_CMP,subtype,sync);
        case 104: return od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_CMP,subtype,sync);
        case 105: return od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_CMP,subtype,sync);
        case 106: return od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_CMP,subtype,sync);
        case 107: return od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_CMP,subtype,sync);
        case 121: return od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_SBS,subtype,sync);
        case 122: return od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_SBS,subtype,sync);
        case 123: return od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_SBS,subtype,sync);
        case 124: return od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_SBS,subtype,sync);
        case 125: return od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_SBS,subtype,sync);
        case 126: return od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_SBS,subtype,sync);
        case 127: return od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_SBS,subtype,sync);
    }
    trace(2,"rtcm3 4076: unsupported message subtype=%d\n",subtype);
    return 0;
}
/* encode RTCM ver.3 message -------------------------------------------------*/
extern int encode_rtcm3(rtcm_t *rtcm, int type, int subtype, int sync)
{
    int ret=0;
    
    trace(3,"encode_rtcm3: type=%d subtype=%d sync=%d\n",type,subtype,sync);
    
    switch (type) {
        case 1001: ret=od_rtk_rtcm3e_encode_type1001(rtcm,sync);     break;
        case 1002: ret=od_rtk_rtcm3e_encode_type1002(rtcm,sync);     break;
        case 1003: ret=od_rtk_rtcm3e_encode_type1003(rtcm,sync);     break;
        case 1004: ret=od_rtk_rtcm3e_encode_type1004(rtcm,sync);     break;
        case 1005: ret=od_rtk_rtcm3e_encode_type1005(rtcm,sync);     break;
        case 1006: ret=od_rtk_rtcm3e_encode_type1006(rtcm,sync);     break;
        case 1007: ret=od_rtk_rtcm3e_encode_type1007(rtcm,sync);     break;
        case 1008: ret=od_rtk_rtcm3e_encode_type1008(rtcm,sync);     break;
        case 1009: ret=od_rtk_rtcm3e_encode_type1009(rtcm,sync);     break;
        case 1010: ret=od_rtk_rtcm3e_encode_type1010(rtcm,sync);     break;
        case 1011: ret=od_rtk_rtcm3e_encode_type1011(rtcm,sync);     break;
        case 1012: ret=od_rtk_rtcm3e_encode_type1012(rtcm,sync);     break;
        case 1019: ret=od_rtk_rtcm3e_encode_type1019(rtcm,sync);     break;
        case 1020: ret=od_rtk_rtcm3e_encode_type1020(rtcm,sync);     break;
        case 1033: ret=od_rtk_rtcm3e_encode_type1033(rtcm,sync);     break;
        case 1041: ret=od_rtk_rtcm3e_encode_type1041(rtcm,sync);     break;
        case 1042: ret=od_rtk_rtcm3e_encode_type1042(rtcm,sync);     break;
        case 1044: ret=od_rtk_rtcm3e_encode_type1044(rtcm,sync);     break;
        case 1045: ret=od_rtk_rtcm3e_encode_type1045(rtcm,sync);     break;
        case 1046: ret=od_rtk_rtcm3e_encode_type1046(rtcm,sync);     break;
        case   63: ret=od_rtk_rtcm3e_encode_type63  (rtcm,sync);     break; /* draft */
        case 1057: ret=od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_GPS,0,sync); break;
        case 1058: ret=od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_GPS,0,sync); break;
        case 1059: ret=od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_GPS,0,sync); break;
        case 1060: ret=od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_GPS,0,sync); break;
        case 1061: ret=od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_GPS,0,sync); break;
        case 1062: ret=od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_GPS,0,sync); break;
        case 1063: ret=od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_GLO,0,sync); break;
        case 1064: ret=od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_GLO,0,sync); break;
        case 1065: ret=od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_GLO,0,sync); break;
        case 1066: ret=od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_GLO,0,sync); break;
        case 1067: ret=od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_GLO,0,sync); break;
        case 1068: ret=od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_GLO,0,sync); break;
        case 1071: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_GPS,sync); break;
        case 1072: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_GPS,sync); break;
        case 1073: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_GPS,sync); break;
        case 1074: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_GPS,sync); break;
        case 1075: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_GPS,sync); break;
        case 1076: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_GPS,sync); break;
        case 1077: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_GPS,sync); break;
        case 1081: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_GLO,sync); break;
        case 1082: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_GLO,sync); break;
        case 1083: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_GLO,sync); break;
        case 1084: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_GLO,sync); break;
        case 1085: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_GLO,sync); break;
        case 1086: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_GLO,sync); break;
        case 1087: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_GLO,sync); break;
        case 1091: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_GAL,sync); break;
        case 1092: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_GAL,sync); break;
        case 1093: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_GAL,sync); break;
        case 1094: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_GAL,sync); break;
        case 1095: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_GAL,sync); break;
        case 1096: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_GAL,sync); break;
        case 1097: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_GAL,sync); break;
        case 1101: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_SBS,sync); break;
        case 1102: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_SBS,sync); break;
        case 1103: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_SBS,sync); break;
        case 1104: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_SBS,sync); break;
        case 1105: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_SBS,sync); break;
        case 1106: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_SBS,sync); break;
        case 1107: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_SBS,sync); break;
        case 1111: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_QZS,sync); break;
        case 1112: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_QZS,sync); break;
        case 1113: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_QZS,sync); break;
        case 1114: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_QZS,sync); break;
        case 1115: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_QZS,sync); break;
        case 1116: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_QZS,sync); break;
        case 1117: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_QZS,sync); break;
        case 1121: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_CMP,sync); break;
        case 1122: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_CMP,sync); break;
        case 1123: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_CMP,sync); break;
        case 1124: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_CMP,sync); break;
        case 1125: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_CMP,sync); break;
        case 1126: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_CMP,sync); break;
        case 1127: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_CMP,sync); break;
        case 1131: ret=od_rtk_rtcm3e_encode_msm1(rtcm,SYS_IRN,sync); break;
        case 1132: ret=od_rtk_rtcm3e_encode_msm2(rtcm,SYS_IRN,sync); break;
        case 1133: ret=od_rtk_rtcm3e_encode_msm3(rtcm,SYS_IRN,sync); break;
        case 1134: ret=od_rtk_rtcm3e_encode_msm4(rtcm,SYS_IRN,sync); break;
        case 1135: ret=od_rtk_rtcm3e_encode_msm5(rtcm,SYS_IRN,sync); break;
        case 1136: ret=od_rtk_rtcm3e_encode_msm6(rtcm,SYS_IRN,sync); break;
        case 1137: ret=od_rtk_rtcm3e_encode_msm7(rtcm,SYS_IRN,sync); break;
        case 1230: ret=od_rtk_rtcm3e_encode_type1230(rtcm,sync);     break;
        case 1240: ret=od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_GAL,0,sync); break; /* draft */
        case 1241: ret=od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_GAL,0,sync); break; /* draft */
        case 1242: ret=od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_GAL,0,sync); break; /* draft */
        case 1243: ret=od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_GAL,0,sync); break; /* draft */
        case 1244: ret=od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_GAL,0,sync); break; /* draft */
        case 1245: ret=od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_GAL,0,sync); break; /* draft */
        case 1246: ret=od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_QZS,0,sync); break; /* draft */
        case 1247: ret=od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_QZS,0,sync); break; /* draft */
        case 1248: ret=od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_QZS,0,sync); break; /* draft */
        case 1249: ret=od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_QZS,0,sync); break; /* draft */
        case 1250: ret=od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_QZS,0,sync); break; /* draft */
        case 1251: ret=od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_QZS,0,sync); break; /* draft */
        case 1252: ret=od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_SBS,0,sync); break; /* draft */
        case 1253: ret=od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_SBS,0,sync); break; /* draft */
        case 1254: ret=od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_SBS,0,sync); break; /* draft */
        case 1255: ret=od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_SBS,0,sync); break; /* draft */
        case 1256: ret=od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_SBS,0,sync); break; /* draft */
        case 1257: ret=od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_SBS,0,sync); break; /* draft */
        case 1258: ret=od_rtk_rtcm3e_encode_ssr1(rtcm,SYS_CMP,0,sync); break; /* draft */
        case 1259: ret=od_rtk_rtcm3e_encode_ssr2(rtcm,SYS_CMP,0,sync); break; /* draft */
        case 1260: ret=od_rtk_rtcm3e_encode_ssr3(rtcm,SYS_CMP,0,sync); break; /* draft */
        case 1261: ret=od_rtk_rtcm3e_encode_ssr4(rtcm,SYS_CMP,0,sync); break; /* draft */
        case 1262: ret=od_rtk_rtcm3e_encode_ssr5(rtcm,SYS_CMP,0,sync); break; /* draft */
        case 1263: ret=od_rtk_rtcm3e_encode_ssr6(rtcm,SYS_CMP,0,sync); break; /* draft */
        case   11: ret=od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_GPS,0,sync); break; /* tentative */
        case   12: ret=od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_GAL,0,sync); break; /* tentative */
        case   13: ret=od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_QZS,0,sync); break; /* tentative */
        case   14: ret=od_rtk_rtcm3e_encode_ssr7(rtcm,SYS_CMP,0,sync); break; /* tentative */
        case 4073: ret=od_rtk_rtcm3e_encode_type4073(rtcm,subtype,sync); break;
        case 4076: ret=od_rtk_rtcm3e_encode_type4076(rtcm,subtype,sync); break;
    }
    if (ret>0) {
        if      (1001<=type&&type<=1299) rtcm->nmsg3[type-1000]++; /*   1-299 */
        else if (4070<=type&&type<=4099) rtcm->nmsg3[type-3770]++; /* 300-329 */
        else rtcm->nmsg3[0]++; /* other */
    }
    return ret;
}


#pragma pop_macro("ROUND_U")
#pragma pop_macro("ROUND")
#pragma pop_macro("RANGE_MS")
#pragma pop_macro("PRUNIT_GPS")
#pragma pop_macro("PRUNIT_GLO")
#pragma pop_macro("P2_66")
#pragma pop_macro("P2_59")
#pragma pop_macro("P2_46")
#pragma pop_macro("P2_41")
#pragma pop_macro("P2_34")
#pragma pop_macro("P2_28")
#pragma pop_macro("P2_10")
#pragma pop_macro("MIN")


/* ===== Embedded ephemeris.c ===== */
#pragma push_macro("COS_5")
#undef COS_5
#pragma push_macro("DEFURASSR")
#undef DEFURASSR
#pragma push_macro("ERREPH_GLO")
#undef ERREPH_GLO
#pragma push_macro("J2_GLO")
#undef J2_GLO
#pragma push_macro("MAXAGESSR")
#undef MAXAGESSR
#pragma push_macro("MAXAGESSR_HRCLK")
#undef MAXAGESSR_HRCLK
#pragma push_macro("MAXCCORSSR")
#undef MAXCCORSSR
#pragma push_macro("MAXECORSSR")
#undef MAXECORSSR
#pragma push_macro("MAX_ITER_KEPLER")
#undef MAX_ITER_KEPLER
#pragma push_macro("MU_CMP")
#undef MU_CMP
#pragma push_macro("MU_GAL")
#undef MU_GAL
#pragma push_macro("MU_GLO")
#undef MU_GLO
#pragma push_macro("MU_GPS")
#undef MU_GPS
#pragma push_macro("OMGE_CMP")
#undef OMGE_CMP
#pragma push_macro("OMGE_GAL")
#undef OMGE_GAL
#pragma push_macro("OMGE_GLO")
#undef OMGE_GLO
#pragma push_macro("RE_GLO")
#undef RE_GLO
#pragma push_macro("RTOL_KEPLER")
#undef RTOL_KEPLER
#pragma push_macro("SIN_5")
#undef SIN_5
#pragma push_macro("SQR")
#undef SQR
#pragma push_macro("STD_BRDCCLK")
#undef STD_BRDCCLK
#pragma push_macro("STD_GAL_NAPA")
#undef STD_GAL_NAPA
#pragma push_macro("TSTEP")
#undef TSTEP

/*------------------------------------------------------------------------------
* ephemeris.c : satellite ephemeris and clock functions
*
*          Copyright (C) 2010-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1] IS-GPS-200K, Navstar GPS Space Segment/Navigation User Interfaces,
*         May 6, 2019
*     [2] Global Navigation Satellite System GLONASS, Interface Control Document
*         Navigational radiosignal In bands L1, L2, (Version 5.1), 2008
*     [3] RTCA/DO-229C, Minimum operational performance standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [4] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*         Code Biases, URA
*     [5] RTCM Paper 012-2009-SC104-528, January 28, 2009 (previous ver of [4])
*     [6] RTCM Paper 012-2009-SC104-582, February 2, 2010 (previous ver of [4])
*     [7] European GNSS (Galileo) Open Service Signal In Space Interface Control
*         Document, Issue 1.3, December, 2016
*     [8] Quasi-Zenith Satellite System Interface Specification Satellite
*         Positioning, Navigation and Timing Service (IS-QZSS-PNT-003), Cabinet
*         Office, November 5, 2018
*     [9] BeiDou navigation satellite system signal in space interface control
*         document open service signal B1I (version 3.0), China Satellite
*         Navigation office, February, 2019
*     [10] RTCM Standard 10403.3, Differential GNSS (Global Navigation
*         Satellite Systems) Services - version 3, October 7, 2016
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.1  moved from rtkcmn.c
*                           added api:
*                               eph2clk(),geph2clk(),seph2clk(),satantoff()
*                               satposs()
*                           changed api:
*                               eph2pos(),geph2pos(),satpos()
*                           deleted api:
*                               satposv(),satposiode()
*           2010/08/26 1.2  add ephemeris option EPHOPT_LEX
*           2010/09/09 1.3  fix problem when precise clock outage
*           2011/01/12 1.4  add api alm2pos()
*                           change api satpos(),satposs()
*                           enable valid unhealthy satellites and output status
*                           fix bug on exception by glonass ephem computation
*           2013/01/10 1.5  support beidou (compass)
*                           use newton's method to solve kepler eq.
*                           update ssr correction algorithm
*           2013/03/20 1.6  fix problem on ssr clock relativitic correction
*           2013/09/01 1.7  support negative pseudorange
*                           fix bug on variance in case of ura ssr = 63
*           2013/11/11 1.8  change constant MAXAGESSR 70.0 -> 90.0
*           2014/10/24 1.9  fix bug on return of var_uraeph() if ura<0||15<ura
*           2014/12/07 1.10 modify MAXDTOE for qzss,gal and bds
*                           test max number of iteration for Kepler
*           2015/08/26 1.11 update RTOL_ELPLER 1E-14 -> 1E-13
*                           set MAX_ITER_KEPLER for alm2pos()
*           2017/04/11 1.12 fix bug on max number of obs data in satposs()
*           2018/10/10 1.13 update reference [7]
*                           support ura value in var_uraeph() for galileo
*                           test eph->flag to recognize beidou geo
*                           add api satseleph() for ephemeris selection
*           2020/11/30 1.14 update references [1],[2],[8],[9] and [10]
*                           add API getseleph()
*                           rename API satseleph() as setseleph()
*                           support NavIC/IRNSS by API satpos() and satposs()
*                           support BDS C59-63 as GEO satellites in eph2pos()
*                           default selection of I/NAV for Galileo ephemeris
*                           no support EPHOPT_LEX by API satpos() and satposs()
*                           unselect Galileo ephemeris with AOD<=0 in seleph()
*                           fix bug on clock iteration in eph2clk(), geph2clk()
*                           fix bug on clock reference time in satpos_ssr()
*                           fix bug on wrong value with ura=15 in var_ura()
*                           use integer types in stdint.h
*-----------------------------------------------------------------------------*/

/* constants and macros ------------------------------------------------------*/

#define SQR(x)   ((x)*(x))

#define RE_GLO   6378136.0        /* radius of earth (m)            ref [2] */
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */

#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */

#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */

#define ERREPH_GLO 5.0            /* error of glonass ephemeris (m) */
#define TSTEP    60.0             /* integration step glonass ephemeris (s) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */

#define DEFURASSR 0.15            /* default accurary of ssr corr (m) */
#define MAXECORSSR 10.0           /* max orbit correction of ssr (m) */
#define MAXCCORSSR (1E-6*CLIGHT)  /* max clock correction of ssr (m) */
#define MAXAGESSR 90.0            /* max age of ssr orbit and clock (s) */
#define MAXAGESSR_HRCLK 10.0      /* max age of ssr high-rate clock (s) */
#define STD_BRDCCLK 30.0          /* error of broadcast clock (m) */
#define STD_GAL_NAPA 500.0        /* error of galileo ephemeris for NAPA (m) */

#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */

/* ephemeris selections ------------------------------------------------------*/
static int od_rtk_ephemeris_eph_sel[]={ /* GPS,GLO,GAL,QZS,BDS,IRN,SBS */
    0,0,0,0,0,0,0
};

/* variance by ura ephemeris -------------------------------------------------*/
static double od_rtk_ephemeris_var_uraeph(int sys, int ura)
{
    const double ura_value[]={   
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0
    };
    if (sys==SYS_GAL) { /* galileo sisa (ref [7] 5.1.11) */
        if (ura<= 49) return SQR(ura*0.01);
        if (ura<= 74) return SQR(0.5+(ura- 50)*0.02);
        if (ura<= 99) return SQR(1.0+(ura- 75)*0.04);
        if (ura<=125) return SQR(2.0+(ura-100)*0.16);
        return SQR(STD_GAL_NAPA);
    }
    else { /* gps ura (ref [1] 20.3.3.3.1.1) */
        return ura<0||14<ura?SQR(6144.0):SQR(ura_value[ura]);
    }
}
/* variance by ura ssr (ref [10] table 3.3-1 DF389) --------------------------*/
static double od_rtk_ephemeris_var_urassr(int ura)
{
    double std;
    if (ura<= 0) return SQR(DEFURASSR);
    if (ura>=63) return SQR(5.4665);
    std=(pow(3.0,(ura>>3)&7)*(1.0+(ura&7)/4.0)-1.0)*1E-3;
    return SQR(std);
}
/* almanac to satellite position and clock bias --------------------------------
* compute satellite position and clock bias with almanac (gps, galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          alm_t *alm       I   almanac
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
* return : none
* notes  : see ref [1],[7],[8]
*-----------------------------------------------------------------------------*/
extern void alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,x,y,sinO,cosO,cosi,mu;
    int n;
    
    trace(4,"alm2pos : time=%s sat=%2d\n",time_str(time,3),alm->sat);
    
    tk=timediff(time,alm->toa);
    
    if (alm->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=0.0;
        return;
    }
    mu=satsys(alm->sat,NULL)==SYS_GAL?MU_GAL:MU_GPS;
    
    M=alm->M0+sqrt(mu/(alm->A*alm->A*alm->A))*tk;
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-alm->e*sin(E)-M)/(1.0-alm->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        trace(2,"alm2pos: kepler iteration overflow sat=%2d\n",alm->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);
    u=atan2(sqrt(1.0-alm->e*alm->e)*sinE,cosE-alm->e)+alm->omg;
    r=alm->A*(1.0-alm->e*cosE);
    i=alm->i0;
    O=alm->OMG0+(alm->OMGd-OMGE)*tk-OMGE*alm->toas;
    x=r*cos(u); y=r*sin(u); sinO=sin(O); cosO=cos(O); cosi=cos(i);
    rs[0]=x*cosO-y*cosi*sinO;
    rs[1]=x*sinO+y*cosi*cosO;
    rs[2]=y*sin(i);
    *dts=alm->f0+alm->f1*tk;
}
/* broadcast ephemeris to satellite clock bias ---------------------------------
* compute satellite clock bias with broadcast ephemeris (gps, galileo, qzss)
* args   : gtime_t time     I   time by satellite clock (gpst)
*          eph_t *eph       I   broadcast ephemeris
* return : satellite clock bias (s) without relativeity correction
* notes  : see ref [1],[7],[8]
*          satellite clock does not include relativity correction and tdg
*-----------------------------------------------------------------------------*/
extern double eph2clk(gtime_t time, const eph_t *eph)
{
    double t,ts;
    int i;
    
    trace(4,"eph2clk : time=%s sat=%2d\n",time_str(time,3),eph->sat);
    
    t=ts=timediff(time,eph->toc);
    
    for (i=0;i<2;i++) {
        t=ts-(eph->f0+eph->f1*t+eph->f2*t*t);
    }
    return eph->f0+eph->f1*t+eph->f2*t*t;
}
/* broadcast ephemeris to satellite position and clock bias --------------------
* compute satellite position and clock bias with broadcast ephemeris (gps,
* galileo, qzss)
* args   : gtime_t time     I   time (gpst)
*          eph_t *eph       I   broadcast ephemeris
*          double *rs       O   satellite position (ecef) {x,y,z} (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [1],[7],[8]
*          satellite clock includes relativity correction without code bias
*          (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern void eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
                    double *var)
{
    double tk,M,E,Ek,sinE,cosE,u,r,i,O,sin2u,cos2u,x,y,sinO,cosO,cosi,mu,omge;
    double xg,yg,zg,sino,coso;
    int n,sys,prn;
    
    trace(4,"eph2pos : time=%s sat=%2d\n",time_str(time,3),eph->sat);
    
    if (eph->A<=0.0) {
        rs[0]=rs[1]=rs[2]=*dts=*var=0.0;
        return;
    }
    tk=timediff(time,eph->toe);
    
    switch ((sys=satsys(eph->sat,&prn))) {
        case SYS_GAL: mu=MU_GAL; omge=OMGE_GAL; break;
        case SYS_CMP: mu=MU_CMP; omge=OMGE_CMP; break;
        default:      mu=MU_GPS; omge=OMGE;     break;
    }
    M=eph->M0+(sqrt(mu/(eph->A*eph->A*eph->A))+eph->deln)*tk;
    
    for (n=0,E=M,Ek=0.0;fabs(E-Ek)>RTOL_KEPLER&&n<MAX_ITER_KEPLER;n++) {
        Ek=E; E-=(E-eph->e*sin(E)-M)/(1.0-eph->e*cos(E));
    }
    if (n>=MAX_ITER_KEPLER) {
        trace(2,"eph2pos: kepler iteration overflow sat=%2d\n",eph->sat);
        return;
    }
    sinE=sin(E); cosE=cos(E);
    
    trace(4,"kepler: sat=%2d e=%8.5f n=%2d del=%10.3e\n",eph->sat,eph->e,n,E-Ek);
    
    u=atan2(sqrt(1.0-eph->e*eph->e)*sinE,cosE-eph->e)+eph->omg;
    r=eph->A*(1.0-eph->e*cosE);
    i=eph->i0+eph->idot*tk;
    sin2u=sin(2.0*u); cos2u=cos(2.0*u);
    u+=eph->cus*sin2u+eph->cuc*cos2u;
    r+=eph->crs*sin2u+eph->crc*cos2u;
    i+=eph->cis*sin2u+eph->cic*cos2u;
    x=r*cos(u); y=r*sin(u); cosi=cos(i);
    
    /* beidou geo satellite */
    if (sys==SYS_CMP&&(prn<=5||prn>=59)) { /* ref [9] table 4-1 */
        O=eph->OMG0+eph->OMGd*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        xg=x*cosO-y*cosi*sinO;
        yg=x*sinO+y*cosi*cosO;
        zg=y*sin(i);
        sino=sin(omge*tk); coso=cos(omge*tk);
        rs[0]= xg*coso+yg*sino*COS_5+zg*sino*SIN_5;
        rs[1]=-xg*sino+yg*coso*COS_5+zg*coso*SIN_5;
        rs[2]=-yg*SIN_5+zg*COS_5;
    }
    else {
        O=eph->OMG0+(eph->OMGd-omge)*tk-omge*eph->toes;
        sinO=sin(O); cosO=cos(O);
        rs[0]=x*cosO-y*cosi*sinO;
        rs[1]=x*sinO+y*cosi*cosO;
        rs[2]=y*sin(i);
    }
    tk=timediff(time,eph->toc);
    *dts=eph->f0+eph->f1*tk+eph->f2*tk*tk;
    
    /* relativity correction */
    *dts-=2.0*sqrt(mu*eph->A)*eph->e*sinE/SQR(CLIGHT);
    
    /* position and clock error variance */
    *var=od_rtk_ephemeris_var_uraeph(sys,eph->sva);
}
/* glonass orbit differential equations --------------------------------------*/
static void od_rtk_ephemeris_deq(const double *x, double *xdot, const double *acc)
{
    double a,b,c,r2=dot(x,x,3),r3=r2*sqrt(r2),omg2=SQR(OMGE_GLO);
    
    if (r2<=0.0) {
        xdot[0]=xdot[1]=xdot[2]=xdot[3]=xdot[4]=xdot[5]=0.0;
        return;
    }
    /* ref [2] A.3.1.2 with bug fix for xdot[4],xdot[5] */
    a=1.5*J2_GLO*MU_GLO*SQR(RE_GLO)/r2/r3; /* 3/2*J2*mu*Ae^2/r^5 */
    b=5.0*x[2]*x[2]/r2;                    /* 5*z^2/r^2 */
    c=-MU_GLO/r3-a*(1.0-b);                /* -mu/r^3-a(1-b) */
    xdot[0]=x[3]; xdot[1]=x[4]; xdot[2]=x[5];
    xdot[3]=(c+omg2)*x[0]+2.0*OMGE_GLO*x[4]+acc[0];
    xdot[4]=(c+omg2)*x[1]-2.0*OMGE_GLO*x[3]+acc[1];
    xdot[5]=(c-2.0*a)*x[2]+acc[2];
}
/* glonass position and velocity by numerical integration --------------------*/
static void od_rtk_ephemeris_glorbit(double t, double *x, const double *acc)
{
    double k1[6],k2[6],k3[6],k4[6],w[6];
    int i;
    
    od_rtk_ephemeris_deq(x,k1,acc); for (i=0;i<6;i++) w[i]=x[i]+k1[i]*t/2.0;
    od_rtk_ephemeris_deq(w,k2,acc); for (i=0;i<6;i++) w[i]=x[i]+k2[i]*t/2.0;
    od_rtk_ephemeris_deq(w,k3,acc); for (i=0;i<6;i++) w[i]=x[i]+k3[i]*t;
    od_rtk_ephemeris_deq(w,k4,acc);
    for (i=0;i<6;i++) x[i]+=(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*t/6.0;
}
/* glonass ephemeris to satellite clock bias -----------------------------------
* compute satellite clock bias with glonass ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          geph_t *geph     I   glonass ephemeris
* return : satellite clock bias (s)
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern double geph2clk(gtime_t time, const geph_t *geph)
{
    double t,ts;
    int i;
    
    trace(4,"geph2clk: time=%s sat=%2d\n",time_str(time,3),geph->sat);
    
    t=ts=timediff(time,geph->toe);
    
    for (i=0;i<2;i++) {
        t=ts-(-geph->taun+geph->gamn*t);
    }
    return -geph->taun+geph->gamn*t;
}
/* glonass ephemeris to satellite position and clock bias ----------------------
* compute satellite position and clock bias with glonass ephemeris
* args   : gtime_t time     I   time (gpst)
*          geph_t *geph     I   glonass ephemeris
*          double *rs       O   satellite position {x,y,z} (ecef) (m)
*          double *dts      O   satellite clock bias (s)
*          double *var      O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [2]
*-----------------------------------------------------------------------------*/
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                     double *var)
{
    double t,tt,x[6];
    int i;
    
    trace(4,"geph2pos: time=%s sat=%2d\n",time_str(time,3),geph->sat);
    
    t=timediff(time,geph->toe);
    
    *dts=-geph->taun+geph->gamn*t;
    
    for (i=0;i<3;i++) {
        x[i  ]=geph->pos[i];
        x[i+3]=geph->vel[i];
    }
    for (tt=t<0.0?-TSTEP:TSTEP;fabs(t)>1E-9;t-=tt) {
        if (fabs(t)<TSTEP) tt=t;
        od_rtk_ephemeris_glorbit(tt,x,geph->acc);
    }
    for (i=0;i<3;i++) rs[i]=x[i];
    
    *var=SQR(ERREPH_GLO);
}
/* sbas ephemeris to satellite clock bias --------------------------------------
* compute satellite clock bias with sbas ephemeris
* args   : gtime_t time     I   time by satellite clock (gpst)
*          seph_t *seph     I   sbas ephemeris
* return : satellite clock bias (s)
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern double seph2clk(gtime_t time, const seph_t *seph)
{
    double t;
    int i;
    
    trace(4,"seph2clk: time=%s sat=%2d\n",time_str(time,3),seph->sat);
    
    t=timediff(time,seph->t0);
    
    for (i=0;i<2;i++) {
        t-=seph->af0+seph->af1*t;
    }
    return seph->af0+seph->af1*t;
}
/* sbas ephemeris to satellite position and clock bias -------------------------
* compute satellite position and clock bias with sbas ephemeris
* args   : gtime_t time     I   time (gpst)
*          seph_t  *seph    I   sbas ephemeris
*          double  *rs      O   satellite position {x,y,z} (ecef) (m)
*          double  *dts     O   satellite clock bias (s)
*          double  *var     O   satellite position and clock variance (m^2)
* return : none
* notes  : see ref [3]
*-----------------------------------------------------------------------------*/
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                     double *var)
{
    double t;
    int i;
    
    trace(4,"seph2pos: time=%s sat=%2d\n",time_str(time,3),seph->sat);
    
    t=timediff(time,seph->t0);
    
    for (i=0;i<3;i++) {
        rs[i]=seph->pos[i]+seph->vel[i]*t+seph->acc[i]*t*t/2.0;
    }
    *dts=seph->af0+seph->af1*t;
    
    *var=od_rtk_ephemeris_var_uraeph(SYS_SBS,seph->sva);
}
/* select ephememeris --------------------------------------------------------*/
static eph_t *od_rtk_ephemeris_seleph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax,tmin;
    int i,j=-1,sys,sel;
    
    trace(4,"seleph  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);
    
    sys=satsys(sat,NULL);
    switch (sys) {
        case SYS_GPS: tmax=MAXDTOE+1.0    ; sel=od_rtk_ephemeris_eph_sel[0]; break;
        case SYS_GAL: tmax=MAXDTOE_GAL    ; sel=od_rtk_ephemeris_eph_sel[2]; break;
        case SYS_QZS: tmax=MAXDTOE_QZS+1.0; sel=od_rtk_ephemeris_eph_sel[3]; break;
        case SYS_CMP: tmax=MAXDTOE_CMP+1.0; sel=od_rtk_ephemeris_eph_sel[4]; break;
        case SYS_IRN: tmax=MAXDTOE_IRN+1.0; sel=od_rtk_ephemeris_eph_sel[5]; break;
        default: tmax=MAXDTOE+1.0; break;
    }
    tmin=tmax+1.0;
    
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        if (iode>=0&&nav->eph[i].iode!=iode) continue;
        if (sys==SYS_GAL) {
            sel=getseleph(SYS_GAL);
            if (sel==0&&!(nav->eph[i].code&(1<<9))) continue; /* I/NAV */
            if (sel==1&&!(nav->eph[i].code&(1<<8))) continue; /* F/NAV */
            if (timediff(nav->eph[i].toe,time)>=0.0) continue; /* AOD<=0 */
        }
        if ((t=fabs(timediff(nav->eph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->eph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        trace(3,"no broadcast ephemeris: %s sat=%2d iode=%3d\n",
              time_str(time,0),sat,iode);
        return NULL;
    }
    return nav->eph+j;
}
/* select glonass ephememeris ------------------------------------------------*/
static geph_t *od_rtk_ephemeris_selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1;
    
    trace(4,"selgeph : time=%s sat=%2d iode=%2d\n",time_str(time,3),sat,iode);
    
    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) {
        trace(3,"no glonass ephemeris  : %s sat=%2d iode=%2d\n",time_str(time,0),
              sat,iode);
        return NULL;
    }
    return nav->geph+j;
}
/* select sbas ephememeris ---------------------------------------------------*/
static seph_t *od_rtk_ephemeris_selseph(gtime_t time, int sat, const nav_t *nav)
{
    double t,tmax=MAXDTOE_SBS,tmin=tmax+1.0;
    int i,j=-1;
    
    trace(4,"selseph : time=%s sat=%2d\n",time_str(time,3),sat);
    
    for (i=0;i<nav->ns;i++) {
        if (nav->seph[i].sat!=sat) continue;
        if ((t=fabs(timediff(nav->seph[i].t0,time)))>tmax) continue;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (j<0) {
        trace(3,"no sbas ephemeris     : %s sat=%2d\n",time_str(time,0),sat);
        return NULL;
    }
    return nav->seph+j;
}
/* satellite clock with broadcast ephemeris ----------------------------------*/
static int od_rtk_ephemeris_ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  double *dts)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    int sys;
    
    trace(4,"ephclk  : time=%s sat=%2d\n",time_str(time,3),sat);
    
    sys=satsys(sat,NULL);
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=od_rtk_ephemeris_seleph(teph,sat,-1,nav))) return 0;
        *dts=eph2clk(time,eph);
    }
    else if (sys==SYS_GLO) {
        if (!(geph=od_rtk_ephemeris_selgeph(teph,sat,-1,nav))) return 0;
        *dts=geph2clk(time,geph);
    }
    else if (sys==SYS_SBS) {
        if (!(seph=od_rtk_ephemeris_selseph(teph,sat,nav))) return 0;
        *dts=seph2clk(time,seph);
    }
    else return 0;
    
    return 1;
}
/* satellite position and clock by broadcast ephemeris -----------------------*/
static int od_rtk_ephemeris_ephpos(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  int iode, double *rs, double *dts, double *var, int *svh)
{
    eph_t  *eph;
    geph_t *geph;
    seph_t *seph;
    double rst[3],dtst[1],tt=1E-3;
    int i,sys;
    
    trace(4,"ephpos  : time=%s sat=%2d iode=%d\n",time_str(time,3),sat,iode);
    
    sys=satsys(sat,NULL);
    
    *svh=-1;
    
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP||sys==SYS_IRN) {
        if (!(eph=od_rtk_ephemeris_seleph(teph,sat,iode,nav))) return 0;
        eph2pos(time,eph,rs,dts,var);
        time=timeadd(time,tt);
        eph2pos(time,eph,rst,dtst,var);
        *svh=eph->svh;
    }
    else if (sys==SYS_GLO) {
        if (!(geph=od_rtk_ephemeris_selgeph(teph,sat,iode,nav))) return 0;
        geph2pos(time,geph,rs,dts,var);
        time=timeadd(time,tt);
        geph2pos(time,geph,rst,dtst,var);
        *svh=geph->svh;
    }
    else if (sys==SYS_SBS) {
        if (!(seph=od_rtk_ephemeris_selseph(teph,sat,nav))) return 0;
        seph2pos(time,seph,rs,dts,var);
        time=timeadd(time,tt);
        seph2pos(time,seph,rst,dtst,var);
        *svh=seph->svh;
    }
    else return 0;
    
    /* satellite velocity and clock drift by differential approx */
    for (i=0;i<3;i++) rs[i+3]=(rst[i]-rs[i])/tt;
    dts[1]=(dtst[0]-dts[0])/tt;
    
    return 1;
}
/* satellite position and clock with sbas correction -------------------------*/
static int od_rtk_ephemeris_satpos_sbas(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                        double *rs, double *dts, double *var, int *svh)
{
    const sbssatp_t *sbs;
    int i;
    
    trace(4,"satpos_sbas: time=%s sat=%2d\n",time_str(time,3),sat);
    
    /* search sbas satellite correciton */
    for (i=0;i<nav->sbssat.nsat;i++) {
        sbs=nav->sbssat.sat+i;
        if (sbs->sat==sat) break;
    }
    if (i>=nav->sbssat.nsat) {
        trace(2,"no sbas correction for orbit: %s sat=%2d\n",time_str(time,0),sat);
        od_rtk_ephemeris_ephpos(time,teph,sat,nav,-1,rs,dts,var,svh);
        *svh=-1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!od_rtk_ephemeris_ephpos(time,teph,sat,nav,sbs->lcorr.iode,rs,dts,var,svh)) return 0;
    
    /* sbas satellite correction (long term and fast) */
    if (sbssatcorr(time,sat,nav,rs,dts,var)) return 1;
    *svh=-1;
    return 0;
}
/* satellite position and clock with ssr correction --------------------------*/
static int od_rtk_ephemeris_satpos_ssr(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                      int opt, double *rs, double *dts, double *var, int *svh)
{
    const ssr_t *ssr;
    eph_t *eph;
    double t1,t2,t3,er[3],ea[3],ec[3],rc[3],deph[3],dclk,dant[3]={0},tk;
    int i,sys;
    
    trace(4,"satpos_ssr: time=%s sat=%2d\n",time_str(time,3),sat);
    
    ssr=nav->ssr+sat-1;
    
    if (!ssr->t0[0].time) {
        trace(2,"no ssr orbit correction: %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (!ssr->t0[1].time) {
        trace(2,"no ssr clock correction: %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* inconsistency between orbit and clock correction */
    if (ssr->iod[0]!=ssr->iod[1]) {
        trace(2,"inconsist ssr correction: %s sat=%2d iod=%d %d\n",
              time_str(time,0),sat,ssr->iod[0],ssr->iod[1]);
        *svh=-1;
        return 0;
    }
    t1=timediff(time,ssr->t0[0]);
    t2=timediff(time,ssr->t0[1]);
    t3=timediff(time,ssr->t0[2]);
    
    /* ssr orbit and clock correction (ref [4]) */
    if (fabs(t1)>MAXAGESSR||fabs(t2)>MAXAGESSR) {
        trace(2,"age of ssr error: %s sat=%2d t=%.0f %.0f\n",time_str(time,0),
              sat,t1,t2);
        *svh=-1;
        return 0;
    }
    if (ssr->udi[0]>=1.0) t1-=ssr->udi[0]/2.0;
    if (ssr->udi[1]>=1.0) t2-=ssr->udi[1]/2.0;
    
    for (i=0;i<3;i++) deph[i]=ssr->deph[i]+ssr->ddeph[i]*t1;
    dclk=ssr->dclk[0]+ssr->dclk[1]*t2+ssr->dclk[2]*t2*t2;
    
    /* ssr highrate clock correction (ref [4]) */
    if (ssr->iod[0]==ssr->iod[2]&&ssr->t0[2].time&&fabs(t3)<MAXAGESSR_HRCLK) {
        dclk+=ssr->hrclk;
    }
    if (norm(deph,3)>MAXECORSSR||fabs(dclk)>MAXCCORSSR) {
        trace(3,"invalid ssr correction: %s deph=%.1f dclk=%.1f\n",
              time_str(time,0),norm(deph,3),dclk);
        *svh=-1;
        return 0;
    }
    /* satellite postion and clock by broadcast ephemeris */
    if (!od_rtk_ephemeris_ephpos(time,teph,sat,nav,ssr->iode,rs,dts,var,svh)) return 0;
    
    /* satellite clock for gps, galileo and qzss */
    sys=satsys(sat,NULL);
    if (sys==SYS_GPS||sys==SYS_GAL||sys==SYS_QZS||sys==SYS_CMP) {
        if (!(eph=od_rtk_ephemeris_seleph(teph,sat,ssr->iode,nav))) return 0;
        
        /* satellite clock by clock parameters */
        tk=timediff(time,eph->toc);
        dts[0]=eph->f0+eph->f1*tk+eph->f2*tk*tk;
        dts[1]=eph->f1+2.0*eph->f2*tk;
        
        /* relativity correction */
        dts[0]-=2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
    }
    /* radial-along-cross directions in ecef */
    if (!normv3(rs+3,ea)) return 0;
    cross3(rs,rs+3,rc);
    if (!normv3(rc,ec)) {
        *svh=-1;
        return 0;
    }
    cross3(ea,ec,er);
    
    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rs,sat,nav,dant);
    }
    for (i=0;i<3;i++) {
        rs[i]+=-(er[i]*deph[0]+ea[i]*deph[1]+ec[i]*deph[2])+dant[i];
    }
    /* t_corr = t_sv - (dts(brdc) + dclk(ssr) / CLIGHT) (ref [10] eq.3.12-7) */
    dts[0]+=dclk/CLIGHT;
    
    /* variance by ssr ura */
    *var=od_rtk_ephemeris_var_urassr(ssr->ura);
    
    trace(5,"satpos_ssr: %s sat=%2d deph=%6.3f %6.3f %6.3f er=%6.3f %6.3f %6.3f dclk=%6.3f var=%6.3f\n",
          time_str(time,2),sat,deph[0],deph[1],deph[2],er[0],er[1],er[2],dclk,*var);
    
    return 1;
}
/* satellite position and clock ------------------------------------------------
* compute satellite position, velocity and clock
* args   : gtime_t time     I   time (gpst)
*          gtime_t teph     I   time to select ephemeris (gpst)
*          int    sat       I   satellite number
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   sat position and velocity (ecef)
*                               {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts      O   sat clock {bias,drift} (s|s/s)
*          double *var      O   sat position and clock error variance (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : status (1:ok,0:error)
* notes  : satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*-----------------------------------------------------------------------------*/
extern int satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                  const nav_t *nav, double *rs, double *dts, double *var,
                  int *svh)
{
    trace(4,"satpos  : time=%s sat=%2d ephopt=%d\n",time_str(time,3),sat,ephopt);
    
    *svh=0;
    
    switch (ephopt) {
        case EPHOPT_BRDC  : return od_rtk_ephemeris_ephpos     (time,teph,sat,nav,-1,rs,dts,var,svh);
        case EPHOPT_SBAS  : return od_rtk_ephemeris_satpos_sbas(time,teph,sat,nav,   rs,dts,var,svh);
        case EPHOPT_SSRAPC: return od_rtk_ephemeris_satpos_ssr (time,teph,sat,nav, 0,rs,dts,var,svh);
        case EPHOPT_SSRCOM: return od_rtk_ephemeris_satpos_ssr (time,teph,sat,nav, 1,rs,dts,var,svh);
        case EPHOPT_PREC  :
            if (!peph2pos(time,sat,nav,1,rs,dts,var)) break; else return 1;
    }
    *svh=-1;
    return 0;
}
/* satellite positions and clocks ----------------------------------------------
* compute satellite positions, velocities and clocks
* args   : gtime_t teph     I   time to select ephemeris (gpst)
*          obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          int    ephopt    I   ephemeris option (EPHOPT_???)
*          double *rs       O   satellite positions and velocities (ecef)
*          double *dts      O   satellite clocks
*          double *var      O   sat position and clock error variances (m^2)
*          int    *svh      O   sat health flag (-1:correction not available)
* return : none
* notes  : rs [(0:2)+i*6]= obs[i] sat position {x,y,z} (m)
*          rs [(3:5)+i*6]= obs[i] sat velocity {vx,vy,vz} (m/s)
*          dts[(0:1)+i*2]= obs[i] sat clock {bias,drift} (s|s/s)
*          var[i]        = obs[i] sat position and clock error variance (m^2)
*          svh[i]        = obs[i] sat health flag
*          if no navigation data, set 0 to rs[], dts[], var[] and svh[]
*          satellite position and clock are values at signal transmission time
*          satellite position is referenced to antenna phase center
*          satellite clock does not include code bias correction (tgd or bgd)
*          any pseudorange and broadcast ephemeris are always needed to get
*          signal transmission time
*-----------------------------------------------------------------------------*/
extern void satposs(gtime_t teph, const obsd_t *obs, int n, const nav_t *nav,
                    int ephopt, double *rs, double *dts, double *var, int *svh)
{
    gtime_t time[2*MAXOBS]={{0}};
    double dt,pr;
    int i,j;
    
    trace(3,"satposs : teph=%s n=%d ephopt=%d\n",time_str(teph,3),n,ephopt);
    
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        for (j=0;j<6;j++) rs [j+i*6]=0.0;
        for (j=0;j<2;j++) dts[j+i*2]=0.0;
        var[i]=0.0; svh[i]=0;
        
        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;
        
        if (j>=NFREQ) {
            trace(3,"no pseudorange %s sat=%2d\n",time_str(obs[i].time,3),obs[i].sat);
            continue;
        }
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);
        
        /* satellite clock bias by broadcast ephemeris */
        if (!od_rtk_ephemeris_ephclk(time[i],teph,obs[i].sat,nav,&dt)) {
            trace(3,"no broadcast clock %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
            continue;
        }
        time[i]=timeadd(time[i],-dt);
        
        /* satellite position and clock at transmission time */
        if (!satpos(time[i],teph,obs[i].sat,ephopt,nav,rs+i*6,dts+i*2,var+i,
                    svh+i)) {
            trace(3,"no ephemeris %s sat=%2d\n",time_str(time[i],3),obs[i].sat);
            continue;
        }
        /* if no precise clock available, use broadcast clock instead */
        if (dts[i*2]==0.0) {
            if (!od_rtk_ephemeris_ephclk(time[i],teph,obs[i].sat,nav,dts+i*2)) continue;
            dts[1+i*2]=0.0;
            *var=SQR(STD_BRDCCLK);
        }
    }
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        trace(4,"%s sat=%2d rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
              time_str(time[i],6),obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],
              dts[i*2]*1E9,var[i],svh[i]);
    }
}
/* set selected satellite ephemeris --------------------------------------------
* Set selected satellite ephemeris for multiple ones like LNAV - CNAV, I/NAV -
* F/NAV. Call it before calling satpos(),satposs() to use unselected one.
* args   : int    sys       I   satellite system (SYS_???)
*          int    sel       I   selection of ephemeris
*                                 GPS,QZS : 0:LNAV ,1:CNAV  (default: LNAV)
*                                 GAL     : 0:I/NAV,1:F/NAV (default: I/NAV)
*                                 others  : undefined
* return : none
* notes  : default ephemeris selection for galileo is any.
*-----------------------------------------------------------------------------*/
extern void setseleph(int sys, int sel)
{
    switch (sys) {
        case SYS_GPS: od_rtk_ephemeris_eph_sel[0]=sel; break;
        case SYS_GLO: od_rtk_ephemeris_eph_sel[1]=sel; break;
        case SYS_GAL: od_rtk_ephemeris_eph_sel[2]=sel; break;
        case SYS_QZS: od_rtk_ephemeris_eph_sel[3]=sel; break;
        case SYS_CMP: od_rtk_ephemeris_eph_sel[4]=sel; break;
        case SYS_IRN: od_rtk_ephemeris_eph_sel[5]=sel; break;
        case SYS_SBS: od_rtk_ephemeris_eph_sel[6]=sel; break;
    }
}
/* get selected satellite ephemeris -------------------------------------------
* Get the selected satellite ephemeris.
* args   : int    sys       I   satellite system (SYS_???)
* return : selected ephemeris
*            refer setseleph()
*-----------------------------------------------------------------------------*/
extern int getseleph(int sys)
{
    switch (sys) {
        case SYS_GPS: return od_rtk_ephemeris_eph_sel[0];
        case SYS_GLO: return od_rtk_ephemeris_eph_sel[1];
        case SYS_GAL: return od_rtk_ephemeris_eph_sel[2];
        case SYS_QZS: return od_rtk_ephemeris_eph_sel[3];
        case SYS_CMP: return od_rtk_ephemeris_eph_sel[4];
        case SYS_IRN: return od_rtk_ephemeris_eph_sel[5];
        case SYS_SBS: return od_rtk_ephemeris_eph_sel[6];
    }
    return 0;
}


#pragma pop_macro("TSTEP")
#pragma pop_macro("STD_GAL_NAPA")
#pragma pop_macro("STD_BRDCCLK")
#pragma pop_macro("SQR")
#pragma pop_macro("SIN_5")
#pragma pop_macro("RTOL_KEPLER")
#pragma pop_macro("RE_GLO")
#pragma pop_macro("OMGE_GLO")
#pragma pop_macro("OMGE_GAL")
#pragma pop_macro("OMGE_CMP")
#pragma pop_macro("MU_GPS")
#pragma pop_macro("MU_GLO")
#pragma pop_macro("MU_GAL")
#pragma pop_macro("MU_CMP")
#pragma pop_macro("MAX_ITER_KEPLER")
#pragma pop_macro("MAXECORSSR")
#pragma pop_macro("MAXCCORSSR")
#pragma pop_macro("MAXAGESSR_HRCLK")
#pragma pop_macro("MAXAGESSR")
#pragma pop_macro("J2_GLO")
#pragma pop_macro("ERREPH_GLO")
#pragma pop_macro("DEFURASSR")
#pragma pop_macro("COS_5")


/* ===== Embedded ionex.c ===== */
#pragma push_macro("MIN_EL")
#undef MIN_EL
#pragma push_macro("MIN_HGT")
#undef MIN_HGT
#pragma push_macro("SQR")
#undef SQR
#pragma push_macro("VAR_NOTEC")
#undef VAR_NOTEC

/*------------------------------------------------------------------------------
* ionex.c : ionex functions
*
*          Copyright (C) 2011-2013 by T.TAKASU, All rights reserved.
*
* references:
*     [1] S.Schear, W.Gurtner and J.Feltens, IONEX: The IONosphere Map EXchange
*         Format Version 1, February 25, 1998
*     [2] S.Schaer, R.Markus, B.Gerhard and A.S.Timon, Daily Global Ionosphere
*         Maps based on GPS Carrier Phase Data Routinely producted by CODE
*         Analysis Center, Proceeding of the IGS Analysis Center Workshop, 1996
*
* version : $Revision:$ $Date:$
* history : 2011/03/29 1.0 new
*           2013/03/05 1.1 change api readtec()
*                          fix problem in case of lat>85deg or lat<-85deg
*           2014/02/22 1.2 fix problem on compiled as C++
*-----------------------------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */

/* get index -----------------------------------------------------------------*/
static int od_rtk_ionex_getindex(double value, const double *range)
{
    if (range[2]==0.0) return 0;
    if (range[1]>0.0&&(value<range[0]||range[1]<value)) return -1;
    if (range[1]<0.0&&(value<range[1]||range[0]<value)) return -1;
    return (int)floor((value-range[0])/range[2]+0.5);
}
/* get number of items -------------------------------------------------------*/
static int od_rtk_ionex_nitem(const double *range)
{
    return od_rtk_ionex_getindex(range[1],range)+1;
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
static int od_rtk_ionex_dataindex(int i, int j, int k, const int *ndata)
{
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}
/* add tec data to navigation data -------------------------------------------*/
static tec_t *od_rtk_ionex_addtec(const double *lats, const double *lons, const double *hgts,
                     double rb, nav_t *nav)
{
    tec_t *p,*nav_tec;
    gtime_t time0={0};
    int i,n,ndata[3];
    
    trace(3,"addtec  :\n");
    
    ndata[0]=od_rtk_ionex_nitem(lats);
    ndata[1]=od_rtk_ionex_nitem(lons);
    ndata[2]=od_rtk_ionex_nitem(hgts);
    if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0) return NULL;
    
    if (nav->nt>=nav->ntmax) {
        nav->ntmax+=256;
        if (!(nav_tec=(tec_t *)realloc(nav->tec,sizeof(tec_t)*nav->ntmax))) {
            trace(1,"readionex malloc error ntmax=%d\n",nav->ntmax);
            free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
            return NULL;
        }
        nav->tec=nav_tec;
    }
    p=nav->tec+nav->nt;
    p->time=time0;
    p->rb=rb;
    for (i=0;i<3;i++) {
        p->ndata[i]=ndata[i];
        p->lats[i]=lats[i];
        p->lons[i]=lons[i];
        p->hgts[i]=hgts[i];
    }
    n=ndata[0]*ndata[1]*ndata[2];
    
    if (!(p->data=(double *)malloc(sizeof(double)*n))||
        !(p->rms =(float  *)malloc(sizeof(float )*n))) {
        return NULL;
    }
    for (i=0;i<n;i++) {
        p->data[i]=0.0;
        p->rms [i]=0.0f;
    }
    nav->nt++;
    return p;
}
/* read ionex dcb aux data ----------------------------------------------------*/
static void od_rtk_ionex_readionexdcb(FILE *fp, double *dcb, double *rms)
{
    int i,sat;
    char buff[1024],id[32],*label;
    
    trace(3,"readionexdcb:\n");
    
    for (i=0;i<MAXSAT;i++) dcb[i]=rms[i]=0.0;
    
    while (fgets(buff,sizeof(buff),fp)) {
        if (strlen(buff)<60) continue;
        label=buff+60;
        
        if (strstr(label,"PRN / BIAS / RMS")==label) {
            
            strncpy(id,buff+3,3); id[3]='\0';
            
            if (!(sat=satid2no(id))) {
                trace(2,"ionex invalid satellite: %s\n",id);
                continue;
            }
            dcb[sat-1]=str2num(buff, 6,10);
            rms[sat-1]=str2num(buff,16,10);
        }
        else if (strstr(label,"END OF AUX DATA")==label) break;
    }
}
/* read ionex header ---------------------------------------------------------*/
static double od_rtk_ionex_readionexh(FILE *fp, double *lats, double *lons, double *hgts,
                         double *rb, double *nexp, double *dcb, double *rms)
{
    double ver=0.0;
    char buff[1024],*label;
    
    trace(3,"readionexh:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (strlen(buff)<60) continue;
        label=buff+60;
        
        if (strstr(label,"IONEX VERSION / TYPE")==label) {
            if (buff[20]=='I') ver=str2num(buff,0,8);
        }
        else if (strstr(label,"BASE RADIUS")==label) {
            *rb=str2num(buff,0,8);
        }
        else if (strstr(label,"HGT1 / HGT2 / DHGT")==label) {
            hgts[0]=str2num(buff, 2,6);
            hgts[1]=str2num(buff, 8,6);
            hgts[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LAT1 / LAT2 / DLAT")==label) {
            lats[0]=str2num(buff, 2,6);
            lats[1]=str2num(buff, 8,6);
            lats[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LON1 / LON2 / DLON")==label) {
            lons[0]=str2num(buff, 2,6);
            lons[1]=str2num(buff, 8,6);
            lons[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"EXPONENT")==label) {
            *nexp=str2num(buff,0,6);
        }
        else if (strstr(label,"START OF AUX DATA")==label&&
                 strstr(buff,"DIFFERENTIAL CODE BIASES")) {
            od_rtk_ionex_readionexdcb(fp,dcb,rms);
        }
        else if (strstr(label,"END OF HEADER")==label) {
            return ver;
        }
    }
    return 0.0;
}
/* read ionex body -----------------------------------------------------------*/
static int od_rtk_ionex_readionexb(FILE *fp, const double *lats, const double *lons,
                      const double *hgts, double rb, double nexp, nav_t *nav)
{
    tec_t *p=NULL;
    gtime_t time={0};
    double lat,lon[3],hgt,x;
    int i,j,k,n,m,index,type=0;
    char buff[1024],*label=buff+60;
    
    trace(3,"readionexb:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (strlen(buff)<60) continue;
        
        if (strstr(label,"START OF TEC MAP")==label) {
            if ((p=od_rtk_ionex_addtec(lats,lons,hgts,rb,nav))) type=1;
        }
        else if (strstr(label,"END OF TEC MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"START OF RMS MAP")==label) {
            type=2;
            p=NULL;
        }
        else if (strstr(label,"END OF RMS MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"EPOCH OF CURRENT MAP")==label) {
            if (str2time(buff,0,36,&time)) {
                trace(2,"ionex epoch invalid: %-36.36s\n",buff);
                continue;
            }
            if (type==2) {
                for (i=nav->nt-1;i>=0;i--) {
                    if (fabs(timediff(time,nav->tec[i].time))>=1.0) continue;
                    p=nav->tec+i;
                    break;
                }
            }
            else if (p) p->time=time;
        }
        else if (strstr(label,"LAT/LON1/LON2/DLON/H")==label&&p) {
            lat   =str2num(buff, 2,6);
            lon[0]=str2num(buff, 8,6);
            lon[1]=str2num(buff,14,6);
            lon[2]=str2num(buff,20,6);
            hgt   =str2num(buff,26,6);
            
            i=od_rtk_ionex_getindex(lat,p->lats);
            k=od_rtk_ionex_getindex(hgt,p->hgts);
            n=od_rtk_ionex_nitem(lon);
            
            for (m=0;m<n;m++) {
                if (m%16==0&&!fgets(buff,sizeof(buff),fp)) break;
                
                j=od_rtk_ionex_getindex(lon[0]+lon[2]*m,p->lons);
                if ((index=od_rtk_ionex_dataindex(i,j,k,p->ndata))<0) continue;
                
                if ((x=str2num(buff,m%16*5,5))==9999.0) continue;
                
                if (type==1) p->data[index]=x*pow(10.0,nexp);
                else p->rms[index]=(float)(x*pow(10.0,nexp));
            }
        }
    }
    return 1;
}
/* combine tec grid data -----------------------------------------------------*/
static void od_rtk_ionex_combtec(nav_t *nav)
{
    tec_t tmp;
    int i,j,n=0;
    
    trace(3,"combtec : nav->nt=%d\n",nav->nt);
    
    for (i=0;i<nav->nt-1;i++) {
        for (j=i+1;j<nav->nt;j++) {
            if (timediff(nav->tec[j].time,nav->tec[i].time)<0.0) {
                tmp=nav->tec[i];
                nav->tec[i]=nav->tec[j];
                nav->tec[j]=tmp;
            }
        }
    }
    for (i=0;i<nav->nt;i++) {
        if (i>0&&timediff(nav->tec[i].time,nav->tec[n-1].time)==0.0) {
            free(nav->tec[n-1].data);
            free(nav->tec[n-1].rms );
            nav->tec[n-1]=nav->tec[i];
            continue;
        }
        nav->tec[n++]=nav->tec[i];
    }
    nav->nt=n;
    
    trace(4,"combtec : nav->nt=%d\n",nav->nt);
}
/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : none
* notes  : see ref [1]
*-----------------------------------------------------------------------------*/
extern void readtec(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    double lats[3]={0},lons[3]={0},hgts[3]={0},rb=0.0,nexp=-1.0;
    double dcb[MAXSAT]={0},rms[MAXSAT]={0};
    int i,n;
    char *efiles[MAXEXFILE];
    
    trace(3,"readtec : file=%s\n",file);
    
    /* clear of tec grid data option */
    if (!opt) {
        free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"ionex file open error %s\n",efiles[i]);
            continue;
        }
        /* read ionex header */
        if (od_rtk_ionex_readionexh(fp,lats,lons,hgts,&rb,&nexp,dcb,rms)<=0.0) {
            trace(2,"ionex file format error %s\n",efiles[i]);
            continue;
        }
        /* read ionex body */
        od_rtk_ionex_readionexb(fp,lats,lons,hgts,rb,nexp,nav);
        
        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    /* combine tec grid data */
    if (nav->nt>0) od_rtk_ionex_combtec(nav);
    
    /* P1-P2 dcb */
    for (i=0;i<MAXSAT;i++) {
        nav->cbias[i][0]=CLIGHT*dcb[i]*1E-9; /* ns->m */
    }
}
/* interpolate tec grid data -------------------------------------------------*/
static int od_rtk_ionex_interptec(const tec_t *tec, int k, const double *posp, double *value,
                     double *rms)
{
    double dlat,dlon,a,b,d[4]={0},r[4]={0};
    int i,j,n,index;
    
    trace(3,"interptec: k=%d posp=%.2f %.2f\n",k,posp[0]*R2D,posp[1]*R2D);
    *value=*rms=0.0;
    
    if (tec->lats[2]==0.0||tec->lons[2]==0.0) return 0;
    
    dlat=posp[0]*R2D-tec->lats[0];
    dlon=posp[1]*R2D-tec->lons[0];
    if (tec->lons[2]>0.0) dlon-=floor( dlon/360)*360.0; /*  0<=dlon<360 */
    else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */
    
    a=dlat/tec->lats[2];
    b=dlon/tec->lons[2];
    i=(int)floor(a); a-=i;
    j=(int)floor(b); b-=j;
    
    /* get gridded tec data */
    for (n=0;n<4;n++) {
        if ((index=od_rtk_ionex_dataindex(i+(n%2),j+(n<2?0:1),k,tec->ndata))<0) continue;
        d[n]=tec->data[index];
        r[n]=tec->rms [index];
    }
    if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) {
        
        /* bilinear interpolation (inside of grid) */
        *value=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
        *rms  =(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
    }
    /* nearest-neighbour extrapolation (outside of grid) */
    else if (a<=0.5&&b<=0.5&&d[0]>0.0) {*value=d[0]; *rms=r[0];}
    else if (a> 0.5&&b<=0.5&&d[1]>0.0) {*value=d[1]; *rms=r[1];}
    else if (a<=0.5&&b> 0.5&&d[2]>0.0) {*value=d[2]; *rms=r[2];}
    else if (a> 0.5&&b> 0.5&&d[3]>0.0) {*value=d[3]; *rms=r[3];}
    else {
        i=0;
        for (n=0;n<4;n++) if (d[n]>0.0) {i++; *value+=d[n]; *rms+=r[n];}
        if(i==0) return 0;
        *value/=i; *rms/=i;
    }
    return 1;
}
/* ionosphere delay by tec grid data -----------------------------------------*/
static int od_rtk_ionex_iondelay(gtime_t time, const tec_t *tec, const double *pos,
                    const double *azel, int opt, double *delay, double *var)
{
    const double fact=40.30E16/FREQ1/FREQ1; /* tecu->L1 iono (m) */
    double fs,posp[3]={0},vtec,rms,hion,rp;
    int i;
    
    trace(3,"iondelay: time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
          pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
    
    *delay=*var=0.0;
    
    for (i=0;i<tec->ndata[2];i++) { /* for a layer */
        
        hion=tec->hgts[0]+tec->hgts[2]*i;
        
        /* ionospheric pierce point position */
        fs=ionppp(pos,azel,tec->rb,hion,posp);
        
        if (opt&2) {
            /* modified single layer mapping function (M-SLM) ref [2] */
            rp=tec->rb/(tec->rb+hion)*sin(0.9782*(PI/2.0-azel[1]));
            fs=1.0/sqrt(1.0-rp*rp);
        }
        if (opt&1) {
            /* earth rotation correction (sun-fixed coordinate) */
            posp[1]+=2.0*PI*timediff(time,tec->time)/86400.0;
        }
        /* interpolate tec grid data */
        if (!od_rtk_ionex_interptec(tec,i,posp,&vtec,&rms)) return 0;
        
        *delay+=fact*fs*vtec;
        *var+=fact*fact*fs*fs*rms*rms;
    }
    trace(4,"iondelay: delay=%7.2f std=%6.2f\n",*delay,sqrt(*var));
    
    return 1;
}
/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
extern int iontec(gtime_t time, const nav_t *nav, const double *pos,
                  const double *azel, int opt, double *delay, double *var)
{
    double dels[2],vars[2],a,tt;
    int i,stat[2];
    
    trace(3,"iontec  : time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
          pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
    
    if (azel[1]<MIN_EL||pos[2]<MIN_HGT) {
        *delay=0.0;
        *var=VAR_NOTEC;
        return 1;
    }
    for (i=0;i<nav->nt;i++) {
        if (timediff(nav->tec[i].time,time)>0.0) break;
    }
    if (i==0||i>=nav->nt) {
        trace(2,"%s: tec grid out of period\n",time_str(time,0));
        return 0;
    }
    if ((tt=timediff(nav->tec[i].time,nav->tec[i-1].time))==0.0) {
        trace(2,"tec grid time interval error\n");
        return 0;
    }
    /* ionospheric delay by tec grid data */
    stat[0]=od_rtk_ionex_iondelay(time,nav->tec+i-1,pos,azel,opt,dels  ,vars  );
    stat[1]=od_rtk_ionex_iondelay(time,nav->tec+i  ,pos,azel,opt,dels+1,vars+1);
    
    if (!stat[0]&&!stat[1]) {
        trace(2,"%s: tec grid out of area pos=%6.2f %7.2f azel=%6.1f %5.1f\n",
              time_str(time,0),pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
        return 0;
    }
    if (stat[0]&&stat[1]) { /* linear interpolation by time */
        a=timediff(time,nav->tec[i-1].time)/tt;
        *delay=dels[0]*(1.0-a)+dels[1]*a;
        *var  =vars[0]*(1.0-a)+vars[1]*a;
    }
    else if (stat[0]) { /* nearest-neighbour extrapolation by time */
        *delay=dels[0];
        *var  =vars[0];
    }
    else {
        *delay=dels[1];
        *var  =vars[1];
    }
    trace(3,"iontec  : delay=%5.2f std=%5.2f\n",*delay,sqrt(*var));
    return 1;
}


#pragma pop_macro("VAR_NOTEC")
#pragma pop_macro("SQR")
#pragma pop_macro("MIN_HGT")
#pragma pop_macro("MIN_EL")


/* ===== Embedded pntpos.c ===== */
#pragma push_macro("ERR_BRDCI")
#undef ERR_BRDCI
#pragma push_macro("ERR_CBIAS")
#undef ERR_CBIAS
#pragma push_macro("ERR_ION")
#undef ERR_ION
#pragma push_macro("ERR_SAAS")
#undef ERR_SAAS
#pragma push_macro("ERR_TROP")
#undef ERR_TROP
#pragma push_macro("MAXITR")
#undef MAXITR
#pragma push_macro("MIN_EL")
#undef MIN_EL
#pragma push_macro("NX")
#undef NX
#pragma push_macro("REL_HUMI")
#undef REL_HUMI
#pragma push_macro("SQR")
#undef SQR

/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*           2020/11/30 1.7  support NavIC/IRNSS in pntpos()
*                           no support IONOOPT_LEX option in ioncorr()
*                           improve handling of TGD correction for each system
*                           use E1-E5b for Galileo dual-freq iono-correction
*                           use API sat2freq() to get carrier frequency
*                           add output of velocity estimation error in estvel()
*-----------------------------------------------------------------------------*/

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#if 0 /* enable GPS-QZS time offset estimation */
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

/* pseudorange measurement error variance ------------------------------------*/
static double od_rtk_pntpos_varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    if (el<MIN_EL) el=MIN_EL;
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}
/* get group delay parameter (m) ---------------------------------------------*/
static double od_rtk_pntpos_gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask -------------------------------------------------------------*/
static int od_rtk_pntpos_snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    if (testsnr(0,0,azel[1],obs->SNR[0]*SNR_UNIT,&opt->snrmask)) {
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        if (testsnr(0,1,azel[1],obs->SNR[1]*SNR_UNIT,&opt->snrmask)) return 0;
    }
    return 1;
}
/* psendorange with code bias correction -------------------------------------*/
static double od_rtk_pntpos_prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var)
{
    double P1,P2,gamma,b1,b2;
    int sat,sys;
    
    sat=obs->sat;
    sys=satsys(sat,NULL);
    P1=obs->P[0];
    P2=obs->P[1];
    *var=0.0;
    
    if (P1==0.0||(opt->ionoopt==IONOOPT_IFLC&&P2==0.0)) return 0.0;
    
    /* P1-C1,P2-C2 DCB correction */
    if (sys==SYS_GPS||sys==SYS_GLO) {
        if (obs->code[0]==CODE_L1C) P1+=nav->cbias[sat-1][1]; /* C1->P1 */
        if (obs->code[1]==CODE_L2C) P2+=nav->cbias[sat-1][2]; /* C2->P2 */
    }
    if (opt->ionoopt==IONOOPT_IFLC) { /* dual-frequency */
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2,G1-G2 */
            gamma=SQR(FREQ1/FREQ2);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GLO) { /* G1-G2 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_GAL) { /* E1-E5b */
            gamma=SQR(FREQ1/FREQ7);
            if (getseleph(SYS_GAL)) { /* F/NAV */
                P2-=od_rtk_pntpos_gettgd(sat,nav,0)-od_rtk_pntpos_gettgd(sat,nav,1); /* BGD_E5aE5b */
            }
            return (P2-gamma*P1)/(1.0-gamma);
        }
        else if (sys==SYS_CMP) { /* B1-B2 */
            gamma=SQR(((obs->code[0]==CODE_L2I)?FREQ1_CMP:FREQ1)/FREQ2_CMP);
            if      (obs->code[0]==CODE_L2I) b1=od_rtk_pntpos_gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=od_rtk_pntpos_gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=od_rtk_pntpos_gettgd(sat,nav,2)+od_rtk_pntpos_gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            b2=od_rtk_pntpos_gettgd(sat,nav,1); /* TGD_B2I/B2bI (m) */
            return ((P2-gamma*P1)-(b2-gamma*b1))/(1.0-gamma);
        }
        else if (sys==SYS_IRN) { /* L5-S */
            gamma=SQR(FREQ5/FREQ9);
            return (P2-gamma*P1)/(1.0-gamma);
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
            b1=od_rtk_pntpos_gettgd(sat,nav,0); /* TGD (m) */
            return P1-b1;
        }
        else if (sys==SYS_GLO) { /* G1 */
            gamma=SQR(FREQ1_GLO/FREQ2_GLO);
            b1=od_rtk_pntpos_gettgd(sat,nav,0); /* -dtaun (m) */
            return P1-b1/(gamma-1.0);
        }
        else if (sys==SYS_GAL) { /* E1 */
            if (getseleph(SYS_GAL)) b1=od_rtk_pntpos_gettgd(sat,nav,0); /* BGD_E1E5a */
            else                    b1=od_rtk_pntpos_gettgd(sat,nav,1); /* BGD_E1E5b */
            return P1-b1;
        }
        else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
            if      (obs->code[0]==CODE_L2I) b1=od_rtk_pntpos_gettgd(sat,nav,0); /* TGD_B1I */
            else if (obs->code[0]==CODE_L1P) b1=od_rtk_pntpos_gettgd(sat,nav,2); /* TGD_B1Cp */
            else b1=od_rtk_pntpos_gettgd(sat,nav,2)+od_rtk_pntpos_gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
            return P1-b1;
        }
        else if (sys==SYS_IRN) { /* L5 */
            gamma=SQR(FREQ9/FREQ5);
            b1=od_rtk_pntpos_gettgd(sat,nav,0); /* TGD (m) */
            return P1-gamma*b1;
        }
    }
    return P1;
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* GPS broadcast ionosphere model */
    if (ionoopt==IONOOPT_BRDC) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* SBAS ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        return sbsioncorr(time,nav,pos,azel,ion,var);
    }
    /* IONEX TEC model */
    if (ionoopt==IONOOPT_TEC) {
        return iontec(time,nav,pos,azel,1,ion,var);
    }
    /* QZSS broadcast ionosphere model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* Saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* SBAS (MOPS) troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
static int od_rtk_pntpos_rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   double *v, double *H, double *var, double *azel, int *vsat,
                   double *resp, int *ns)
{
    gtime_t time;
    double r,freq,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P;
    int i,j,nv=0,sat,sys,mask[NX-3]={0};
    
    trace(3,"resprng : n=%d\n",n);
    
    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            trace(2,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
            i++;
            continue;
        }
        /* excluded satellite? */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* geometric distance */
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        
        if (iter>0) {
            /* test elevation mask */
            if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
            
            /* test SNR mask */
            if (!od_rtk_pntpos_snrmask(obs+i,azel+i*2,opt)) continue;
            
            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                continue;
            }
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            dion*=SQR(FREQ1/freq);
            vion*=SQR(FREQ1/freq);
            
            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }
        }
        /* psendorange with code bias correction */
        if ((P=od_rtk_pntpos_prange(obs+i,nav,opt,&vmeas))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        /* design matrix */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }
        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
#if 0 /* enable QZS-GPS time offset estimation */
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* variance of pseudorange error */
        var[nv++]=od_rtk_pntpos_varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;
        
        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int od_rtk_pntpos_valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* Chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        sprintf(msg,"chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        return 0;
    }
    /* large GDOP check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int od_rtk_pntpos_estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                  double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    
    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals (m) */
        nv=od_rtk_pntpos_rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                   &ns);
        
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* weighted by Std */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }
        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            
            /* validate solution */
            if ((stat=od_rtk_pntpos_valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);
    return 0;
}
/* RAIM FDE (failure detection and exclution) -------------------------------*/
static int od_rtk_pntpos_raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                    double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!od_rtk_pntpos_estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* range rate residuals ------------------------------------------------------*/
static int od_rtk_pntpos_resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double err, double *v,
                  double *H)
{
    double freq,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,sig;
    int i,j,nv=0;
    
    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);
        
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e);
        
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* Std of range rate error (m/s) */
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        
        /* range rate residual (m/s) */
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        
        /* design matrix */
        for (j=0;j<4;j++) {
            H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
        }
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void od_rtk_pntpos_estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    double err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,nv;
    
    trace(3,"estvel  : n=%d\n",n);
    
    v=mat(n,1); H=mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
        /* range rate residuals (m/s) */
        if ((nv=od_rtk_pntpos_resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time;
    msg[0]='\0';
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position with pseudorange */
    stat=od_rtk_pntpos_estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=od_rtk_pntpos_raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    }
    /* estimate receiver velocity with Doppler */
    if (stat) {
        od_rtk_pntpos_estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
            ssat[i].snr[0]=0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            ssat[obs[i].sat-1].snr[0]=obs[i].SNR[0];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}


#pragma pop_macro("SQR")
#pragma pop_macro("REL_HUMI")
#pragma pop_macro("NX")
#pragma pop_macro("MIN_EL")
#pragma pop_macro("MAXITR")
#pragma pop_macro("ERR_TROP")
#pragma pop_macro("ERR_SAAS")
#pragma pop_macro("ERR_ION")
#pragma pop_macro("ERR_CBIAS")
#pragma pop_macro("ERR_BRDCI")


/* ===== Embedded preceph.c ===== */
#pragma push_macro("EXTERR_CLK")
#undef EXTERR_CLK
#pragma push_macro("EXTERR_EPH")
#undef EXTERR_EPH
#pragma push_macro("MAXDTE")
#undef MAXDTE
#pragma push_macro("NMAX")
#undef NMAX
#pragma push_macro("SQR")
#undef SQR

/*------------------------------------------------------------------------------
* preceph.c : precise ephemeris and clock functions
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* references :
*     [1] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-c),
*         12 February, 2007
*     [2] J.Ray, W.Gurtner, RINEX Extensions to Handle Clock Information,
*         27 August, 1998
*     [3] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*     [4] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*         Space Technology Library, 2004
*     [5] S.Hilla, The Extended Standard Product 3 Orbit Format (SP3-d),
*         February 21, 2016
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2009/01/18 1.0  new
*           2009/01/31 1.1  fix bug on numerical error to read sp3a ephemeris
*           2009/05/15 1.2  support glonass,galileo,qzs
*           2009/12/11 1.3  support wild-card expansion of file path
*           2010/07/21 1.4  added api:
*                               eci2ecef(),sunmoonpos(),peph2pos(),satantoff(),
*                               readdcb()
*                           changed api:
*                               readsp3()
*                           deleted api:
*                               eph2posp()
*           2010/09/09 1.5  fix problem when precise clock outage
*           2011/01/23 1.6  support qzss satellite code
*           2011/09/12 1.7  fix problem on precise clock outage
*                           move sunmmonpos() to rtkcmn.c
*           2011/12/01 1.8  modify api readsp3()
*                           precede later ephemeris if ephemeris is NULL
*                           move eci2ecef() to rtkcmn.c
*           2013/05/08 1.9  fix bug on computing std-dev of precise clocks
*           2013/11/20 1.10 modify option for api readsp3()
*           2014/04/03 1.11 accept extenstion including sp3,eph,SP3,EPH
*           2014/05/23 1.12 add function to read sp3 velocity records
*                           change api: satantoff()
*           2014/08/31 1.13 add member cov and vco in peph_t sturct
*           2014/10/13 1.14 fix bug on clock error variance in peph2pos()
*           2015/05/10 1.15 add api readfcb()
*                           modify api readdcb()
*           2017/04/11 1.16 fix bug on antenna offset correction in peph2pos()
*           2020/11/30 1.17 support SP3-d [5] to accept more than 85 satellites
*                           support NavIC/IRNSS in API peph2pos()
*                           LC defined GPS/QZS L1-L2, GLO G1-G2, GAL E1-E5b,
*                            BDS B1I-B2I and IRN L5-S for API satantoff()
*                           fix bug on reading SP3 file extension
*-----------------------------------------------------------------------------*/

#define SQR(x)      ((x)*(x))

#define NMAX        10              /* order of polynomial interpolation */
#define MAXDTE      900.0           /* max time difference to ephem time (s) */
#define EXTERR_CLK  1E-3            /* extrapolation error for clock (m/s) */
#define EXTERR_EPH  5E-7            /* extrapolation error for ephem (m/s^2) */

/* satellite code to satellite system ----------------------------------------*/
static int od_rtk_preceph_code2sys(char code)
{
    if (code=='G'||code==' ') return SYS_GPS;
    if (code=='R') return SYS_GLO;
    if (code=='E') return SYS_GAL; /* SP3-d */
    if (code=='J') return SYS_QZS; /* SP3-d */
    if (code=='C') return SYS_CMP; /* SP3-d */
    if (code=='I') return SYS_IRN; /* SP3-d */
    if (code=='L') return SYS_LEO; /* SP3-d */
    return SYS_NONE;
}
/* read SP3 header -----------------------------------------------------------*/
static int od_rtk_preceph_readsp3h(FILE *fp, gtime_t *time, char *type, int *sats,
                    double *bfact, char *tsys)
{
    int i,j,k=0,ns=0,sys,prn;
    char buff[1024];
    
    trace(3,"readsp3h:\n");
    
    for (i=0;;i++) {
        if (!fgets(buff,sizeof(buff),fp)) break;
        
        if (i==0) {
            *type=buff[2];
            if (str2time(buff,3,28,time)) return 0;
        }
        else if (!strncmp(buff,"+ ",2)) { /* satellite id */
            if (ns==0) {
                ns=(int)str2num(buff,4,2);
            }
            for (j=0;j<17&&k<ns;j++) {
                sys=od_rtk_preceph_code2sys(buff[9+3*j]);
                prn=(int)str2num(buff,10+3*j,2);
                if (k<MAXSAT) sats[k++]=satno(sys,prn);
            }
        }
        else if (!strncmp(buff,"++",2)) { /* orbit accuracy */
            continue;
        }
        else if (!strncmp(buff,"%c",2)) { /* time system */
            strncpy(tsys,buff+9,3); tsys[3]='\0';
        }
        else if (!strncmp(buff,"%f",2)&&bfact[0]==0.0) { /* fp base number */
            bfact[0]=str2num(buff, 3,10);
            bfact[1]=str2num(buff,14,12);
        }
        else if (!strncmp(buff,"%i",2)) {
            continue;
        }
        else if (!strncmp(buff,"/*",2)) { /* comment */
            continue;
        }
        else if (!strncmp(buff,"* ",2)) { /* first record */
            /* roll back file pointer */
            fseek(fp,-(long)strlen(buff),SEEK_CUR);
            break;
        }
    }
    return ns;
}
/* add precise ephemeris -----------------------------------------------------*/
static int od_rtk_preceph_addpeph(nav_t *nav, peph_t *peph)
{
    peph_t *nav_peph;
    
    if (nav->ne>=nav->nemax) {
        nav->nemax+=256;
        if (!(nav_peph=(peph_t *)realloc(nav->peph,sizeof(peph_t)*nav->nemax))) {
            trace(1,"readsp3b malloc error n=%d\n",nav->nemax);
            free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
            return 0;
        }
        nav->peph=nav_peph;
    }
    nav->peph[nav->ne++]=*peph;
    return 1;
}
/* read SP3 body -------------------------------------------------------------*/
static void od_rtk_preceph_readsp3b(FILE *fp, char type, int *sats, int ns, double *bfact,
                     char *tsys, int index, int opt, nav_t *nav)
{
    peph_t peph;
    gtime_t time;
    double val,std,base;
    int i,j,sat,sys,prn,n=ns*(type=='P'?1:2),pred_o,pred_c,v;
    char buff[1024];
    
    trace(3,"readsp3b: type=%c ns=%d index=%d opt=%d\n",type,ns,index,opt);
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (!strncmp(buff,"EOF",3)) break;
        
        if (buff[0]!='*'||str2time(buff,3,28,&time)) {
            trace(2,"sp3 invalid epoch %31.31s\n",buff);
            continue;
        }
        if (!strcmp(tsys,"UTC")) time=utc2gpst(time); /* utc->gpst */
        peph.time =time;
        peph.index=index;
        
        for (i=0;i<MAXSAT;i++) {
            for (j=0;j<4;j++) {
                peph.pos[i][j]=0.0;
                peph.std[i][j]=0.0f;
                peph.vel[i][j]=0.0;
                peph.vst[i][j]=0.0f;
            }
            for (j=0;j<3;j++) {
                peph.cov[i][j]=0.0f;
                peph.vco[i][j]=0.0f;
            }
        }
        for (i=pred_o=pred_c=v=0;i<n&&fgets(buff,sizeof(buff),fp);i++) {
            
            if (strlen(buff)<4||(buff[0]!='P'&&buff[0]!='V')) continue;
            
            sys=buff[1]==' '?SYS_GPS:od_rtk_preceph_code2sys(buff[1]);
            prn=(int)str2num(buff,2,2);
            if      (sys==SYS_SBS) prn+=100;
            else if (sys==SYS_QZS) prn+=192; /* extension to sp3-c */
            
            if (!(sat=satno(sys,prn))) continue;
            
            if (buff[0]=='P') {
                pred_c=strlen(buff)>=76&&buff[75]=='P';
                pred_o=strlen(buff)>=80&&buff[79]=='P';
            }
            for (j=0;j<4;j++) {
                
                /* read option for predicted value */
                if (j< 3&&(opt&1)&& pred_o) continue;
                if (j< 3&&(opt&2)&&!pred_o) continue;
                if (j==3&&(opt&1)&& pred_c) continue;
                if (j==3&&(opt&2)&&!pred_c) continue;
                
                val=str2num(buff, 4+j*14,14);
                std=str2num(buff,61+j* 3,j<3?2:3);
                
                if (buff[0]=='P') { /* position */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.pos[sat-1][j]=val*(j<3?1000.0:1E-6);
                        v=1; /* valid epoch */
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.std[sat-1][j]=(float)(pow(base,std)*(j<3?1E-3:1E-12));
                    }
                }
                else if (v) { /* velocity */
                    if (val!=0.0&&fabs(val-999999.999999)>=1E-6) {
                        peph.vel[sat-1][j]=val*(j<3?0.1:1E-10);
                    }
                    if ((base=bfact[j<3?0:1])>0.0&&std>0.0) {
                        peph.vst[sat-1][j]=(float)(pow(base,std)*(j<3?1E-7:1E-16));
                    }
                }
            }
        }
        if (v) {
            if (!od_rtk_preceph_addpeph(nav,&peph)) return;
        }
    }
}
/* compare precise ephemeris -------------------------------------------------*/
static int od_rtk_preceph_cmppeph(const void *p1, const void *p2)
{
    peph_t *q1=(peph_t *)p1,*q2=(peph_t *)p2;
    double tt=timediff(q1->time,q2->time);
    return tt<-1E-9?-1:(tt>1E-9?1:q1->index-q2->index);
}
/* combine precise ephemeris -------------------------------------------------*/
static void od_rtk_preceph_combpeph(nav_t *nav, int opt)
{
    int i,j,k,m;
    
    trace(3,"combpeph: ne=%d\n",nav->ne);
    
    qsort(nav->peph,nav->ne,sizeof(peph_t),od_rtk_preceph_cmppeph);
    
    if (opt&4) return;
    
    for (i=0,j=1;j<nav->ne;j++) {
        
        if (fabs(timediff(nav->peph[i].time,nav->peph[j].time))<1E-9) {
            
            for (k=0;k<MAXSAT;k++) {
                if (norm(nav->peph[j].pos[k],4)<=0.0) continue;
                for (m=0;m<4;m++) nav->peph[i].pos[k][m]=nav->peph[j].pos[k][m];
                for (m=0;m<4;m++) nav->peph[i].std[k][m]=nav->peph[j].std[k][m];
                for (m=0;m<4;m++) nav->peph[i].vel[k][m]=nav->peph[j].vel[k][m];
                for (m=0;m<4;m++) nav->peph[i].vst[k][m]=nav->peph[j].vst[k][m];
            }
        }
        else if (++i<j) nav->peph[i]=nav->peph[j];
    }
    nav->ne=i+1;
    
    trace(4,"combpeph: ne=%d\n",nav->ne);
}
/* read sp3 precise ephemeris file ---------------------------------------------
* read sp3 precise ephemeris/clock files and set them to navigation data
* args   : char   *file       I   sp3-c precise ephemeris file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*          int    opt         I   options (1: only observed + 2: only predicted +
*                                 4: not combined)
* return : none
* notes  : see ref [1]
*          precise ephemeris is appended and combined
*          nav->peph and nav->ne must by properly initialized before calling the
*          function
*          only files with extensions of .sp3, .SP3, .eph* and .EPH* are read
*-----------------------------------------------------------------------------*/
extern void readsp3(const char *file, nav_t *nav, int opt)
{
    FILE *fp;
    gtime_t time={0};
    double bfact[2]={0};
    int i,j,n,ns,sats[MAXSAT]={0};
    char *efiles[MAXEXFILE],*ext,type=' ',tsys[4]="";
    
    trace(3,"readpephs: file=%s\n",file);
    
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=j=0;i<n;i++) {
        if (!(ext=strrchr(efiles[i],'.'))) continue;
        
        if (!strstr(ext,".sp3")&&!strstr(ext,".SP3")&&
            !strstr(ext,".eph")&&!strstr(ext,".EPH")) continue;
        
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"sp3 file open error %s\n",efiles[i]);
            continue;
        }
        /* read sp3 header */
        ns=od_rtk_preceph_readsp3h(fp,&time,&type,sats,bfact,tsys);
        
        /* read sp3 body */
        od_rtk_preceph_readsp3b(fp,type,sats,ns,bfact,tsys,j++,opt,nav);
        
        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    /* combine precise ephemeris */
    if (nav->ne>0) od_rtk_preceph_combpeph(nav,opt);
}
/* read satellite antenna parameters -------------------------------------------
* read satellite antenna parameters
* args   : char   *file       I   antenna parameter file
*          gtime_t time       I   time
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : only support antex format for the antenna parameter file
*-----------------------------------------------------------------------------*/
extern int readsap(const char *file, gtime_t time, nav_t *nav)
{
    pcvs_t pcvs={0};
    pcv_t pcv0={0},*pcv;
    int i;
    
    trace(3,"readsap : file=%s time=%s\n",file,time_str(time,0));
    
    if (!readpcv(file,&pcvs)) return 0;
    
    for (i=0;i<MAXSAT;i++) {
        pcv=searchpcv(i+1,"",time,&pcvs);
        nav->pcvs[i]=pcv?*pcv:pcv0;
    }
    free(pcvs.pcv);
    return 1;
}
/* read DCB parameters file --------------------------------------------------*/
static int od_rtk_preceph_readdcbf(const char *file, nav_t *nav, const sta_t *sta)
{
    FILE *fp;
    double cbias;
    char buff[256],str1[32],str2[32]="";
    int i,j,sat,type=0;
    
    trace(3,"readdcbf: file=%s\n",file);
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"dcb parameters file open error: %s\n",file);
        return 0;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        
        if      (strstr(buff,"DIFFERENTIAL (P1-P2) CODE BIASES")) type=1;
        else if (strstr(buff,"DIFFERENTIAL (P1-C1) CODE BIASES")) type=2;
        else if (strstr(buff,"DIFFERENTIAL (P2-C2) CODE BIASES")) type=3;
        
        if (!type||sscanf(buff,"%s %s",str1,str2)<1) continue;
        
        if ((cbias=str2num(buff,26,9))==0.0) continue;
        
        if (sta&&(!strcmp(str1,"G")||!strcmp(str1,"R"))) { /* receiver DCB */
            for (i=0;i<MAXRCV;i++) {
                if (!strcmp(sta[i].name,str2)) break;
            }
            if (i<MAXRCV) {
                j=!strcmp(str1,"G")?0:1;
                nav->rbias[i][j][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
            }
        }
        else if ((sat=satid2no(str1))) { /* satellite dcb */
            nav->cbias[sat-1][type-1]=cbias*1E-9*CLIGHT; /* ns -> m */
        }
    }
    fclose(fp);
    
    return 1;
}
/* read DCB parameters ---------------------------------------------------------
* read differential code bias (DCB) parameters
* args   : char   *file       I   DCB parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
*          sta_t  *sta        I   station info data to inport receiver DCB
*                                 (NULL: no use)
* return : status (1:ok,0:error)
* notes  : currently only support P1-P2, P1-C1, P2-C2, bias in DCB file
*-----------------------------------------------------------------------------*/
extern int readdcb(const char *file, nav_t *nav, const sta_t *sta)
{
    int i,j,n;
    char *efiles[MAXEXFILE]={0};
    
    trace(3,"readdcb : file=%s\n",file);
    
    for (i=0;i<MAXSAT;i++) for (j=0;j<3;j++) {
        nav->cbias[i][j]=0.0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
        od_rtk_preceph_readdcbf(efiles[i],nav,sta);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    return 1;
}
/* polynomial interpolation by Neville's algorithm ---------------------------*/
static double od_rtk_preceph_interppol(const double *x, double *y, int n)
{
    int i,j;
    
    for (j=1;j<n;j++) {
        for (i=0;i<n-j;i++) {
            y[i]=(x[i+j]*y[i]-x[i]*y[i+1])/(x[i+j]-x[i]);
        }
    }
    return y[0];
}
/* satellite position by precise ephemeris -----------------------------------*/
static int od_rtk_preceph_pephpos(gtime_t time, int sat, const nav_t *nav, double *rs,
                   double *dts, double *vare, double *varc)
{
    double t[NMAX+1],p[3][NMAX+1],c[2],*pos,std=0.0,s[3],sinl,cosl;
    int i,j,k,index;
    
    trace(4,"pephpos : time=%s sat=%2d\n",time_str(time,3),sat);
    
    rs[0]=rs[1]=rs[2]=dts[0]=0.0;
    
    if (nav->ne<NMAX+1||
        timediff(time,nav->peph[0].time)<-MAXDTE||
        timediff(time,nav->peph[nav->ne-1].time)>MAXDTE) {
        trace(3,"no prec ephem %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    /* binary search */
    for (i=0,j=nav->ne-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->peph[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* polynomial interpolation for orbit */
    i=index-(NMAX+1)/2;
    if (i<0) i=0; else if (i+NMAX>=nav->ne) i=nav->ne-NMAX-1;
    
    for (j=0;j<=NMAX;j++) {
        t[j]=timediff(nav->peph[i+j].time,time);
        if (norm(nav->peph[i+j].pos[sat-1],3)<=0.0) {
            trace(3,"prec ephem outage %s sat=%2d\n",time_str(time,0),sat);
            return 0;
        }
    }
    for (j=0;j<=NMAX;j++) {
        pos=nav->peph[i+j].pos[sat-1];
        /* correciton for earh rotation ver.2.4.0 */
        sinl=sin(OMGE*t[j]);
        cosl=cos(OMGE*t[j]);
        p[0][j]=cosl*pos[0]-sinl*pos[1];
        p[1][j]=sinl*pos[0]+cosl*pos[1];
        p[2][j]=pos[2];
    }
    for (i=0;i<3;i++) {
        rs[i]=od_rtk_preceph_interppol(t,p[i],NMAX+1);
    }
    if (vare) {
        for (i=0;i<3;i++) s[i]=nav->peph[index].std[sat-1][i];
        std=norm(s,3);
        
        /* extrapolation error for orbit */
        if      (t[0   ]>0.0) std+=EXTERR_EPH*SQR(t[0   ])/2.0;
        else if (t[NMAX]<0.0) std+=EXTERR_EPH*SQR(t[NMAX])/2.0;
        *vare=SQR(std);
    }
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->peph[index  ].time);
    t[1]=timediff(time,nav->peph[index+1].time);
    c[0]=nav->peph[index  ].pos[sat-1][3];
    c[1]=nav->peph[index+1].pos[sat-1][3];
    
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])!=0.0) {
            std=nav->peph[index].std[sat-1][3]*CLIGHT-EXTERR_CLK*t[0];
        }
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])!=0.0) {
            std=nav->peph[index+1].std[sat-1][3]*CLIGHT+EXTERR_CLK*t[1];
        }
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->peph[index+i].std[sat-1][3]+EXTERR_CLK*fabs(t[i]);
    }
    else {
        dts[0]=0.0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite clock by precise clock ------------------------------------------*/
static int od_rtk_preceph_pephclk(gtime_t time, int sat, const nav_t *nav, double *dts,
                   double *varc)
{
    double t[2],c[2],std;
    int i,j,k,index;
    
    trace(4,"pephclk : time=%s sat=%2d\n",time_str(time,3),sat);
    
    if (nav->nc<2||
        timediff(time,nav->pclk[0].time)<-MAXDTE||
        timediff(time,nav->pclk[nav->nc-1].time)>MAXDTE) {
        trace(3,"no prec clock %s sat=%2d\n",time_str(time,0),sat);
        return 1;
    }
    /* binary search */
    for (i=0,j=nav->nc-1;i<j;) {
        k=(i+j)/2;
        if (timediff(nav->pclk[k].time,time)<0.0) i=k+1; else j=k;
    }
    index=i<=0?0:i-1;
    
    /* linear interpolation for clock */
    t[0]=timediff(time,nav->pclk[index  ].time);
    t[1]=timediff(time,nav->pclk[index+1].time);
    c[0]=nav->pclk[index  ].clk[sat-1][0];
    c[1]=nav->pclk[index+1].clk[sat-1][0];
    
    if (t[0]<=0.0) {
        if ((dts[0]=c[0])==0.0) return 0;
        std=nav->pclk[index].std[sat-1][0]*CLIGHT-EXTERR_CLK*t[0];
    }
    else if (t[1]>=0.0) {
        if ((dts[0]=c[1])==0.0) return 0;
        std=nav->pclk[index+1].std[sat-1][0]*CLIGHT+EXTERR_CLK*t[1];
    }
    else if (c[0]!=0.0&&c[1]!=0.0) {
        dts[0]=(c[1]*t[0]-c[0]*t[1])/(t[0]-t[1]);
        i=t[0]<-t[1]?0:1;
        std=nav->pclk[index+i].std[sat-1][0]*CLIGHT+EXTERR_CLK*fabs(t[i]);
    }
    else {
        trace(3,"prec clock outage %s sat=%2d\n",time_str(time,0),sat);
        return 0;
    }
    if (varc) *varc=SQR(std);
    return 1;
}
/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
* notes  : iono-free LC frequencies defined as follows:
*            GPS/QZSS : L1-L2
*            GLONASS  : G1-G2
*            Galileo  : E1-E5b
*            BDS      : B1I-B2I
*            NavIC    : L5-S
*-----------------------------------------------------------------------------*/
extern void satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
                      double *dant)
{
    const pcv_t *pcv=nav->pcvs+sat-1;
    double ex[3],ey[3],ez[3],es[3],r[3],rsun[3],gmst,erpv[5]={0},freq[2];
    double C1,C2,dant1,dant2;
    int i,sys;
    
    trace(4,"satantoff: time=%s sat=%2d\n",time_str(time,3),sat);
    
    dant[0]=dant[1]=dant[2]=0.0;
    
    /* sun position in ecef */
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,&gmst);
    
    /* unit vectors of satellite fixed coordinates */
    for (i=0;i<3;i++) r[i]=-rs[i];
    if (!normv3(r,ez)) return;
    for (i=0;i<3;i++) r[i]=rsun[i]-rs[i];
    if (!normv3(r,es)) return;
    cross3(ez,es,r);
    if (!normv3(r,ey)) return;
    cross3(ey,ez,ex);
    
    /* iono-free LC coefficients */
    sys=satsys(sat,NULL);
    if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 */
        freq[0]=FREQ1;
        freq[1]=FREQ2;
    }
    else if (sys==SYS_GLO) { /* G1-G2 */
        freq[0]=sat2freq(sat,CODE_L1C,nav);
        freq[1]=sat2freq(sat,CODE_L2C,nav);
    }
    else if (sys==SYS_GAL) { /* E1-E5b */
        freq[0]=FREQ1;
        freq[1]=FREQ7;
    }
    else if (sys==SYS_CMP) { /* B1I-B2I */
        freq[0]=FREQ1_CMP;
        freq[1]=FREQ2_CMP;
    }
    else if (sys==SYS_IRN) { /* B1I-B2I */
        freq[0]=FREQ5;
        freq[1]=FREQ9;
    }
    else return;
    
    C1= SQR(freq[0])/(SQR(freq[0])-SQR(freq[1]));
    C2=-SQR(freq[1])/(SQR(freq[0])-SQR(freq[1]));
    
    /* iono-free LC */
    for (i=0;i<3;i++) {
        dant1=pcv->off[0][0]*ex[i]+pcv->off[0][1]*ey[i]+pcv->off[0][2]*ez[i];
        dant2=pcv->off[1][0]*ex[i]+pcv->off[1][1]*ey[i]+pcv->off[1][2]*ez[i];
        dant[i]=C1*dant1+C2*dant2;
    }
}
/* satellite position/clock by precise ephemeris/clock -------------------------
* compute satellite position/clock with precise ephemeris/clock
* args   : gtime_t time       I   time (gpst)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          int    opt         I   sat postion option
*                                 (0: center of mass, 1: antenna phase center)
*          double *rs         O   sat position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          double *dts        O   sat clock {bias,drift} (s|s/s)
*          double *var        IO  sat position and clock error variance (m)
*                                 (NULL: no output)
* return : status (1:ok,0:error or data outage)
* notes  : clock includes relativistic correction but does not contain code bias
*          before calling the function, nav->peph, nav->ne, nav->pclk and
*          nav->nc must be set by calling readsp3(), readrnx() or readrnxt()
*          if precise clocks are not set, clocks in sp3 are used instead
*-----------------------------------------------------------------------------*/
extern int peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                    double *rs, double *dts, double *var)
{
    gtime_t time_tt;
    double rss[3],rst[3],dtss[1],dtst[1],dant[3]={0},vare=0.0,varc=0.0,tt=1E-3;
    int i;
    
    trace(4,"peph2pos: time=%s sat=%2d opt=%d\n",time_str(time,3),sat,opt);
    
    if (sat<=0||MAXSAT<sat) return 0;
    
    /* satellite position and clock bias */
    if (!od_rtk_preceph_pephpos(time,sat,nav,rss,dtss,&vare,&varc)||
        !od_rtk_preceph_pephclk(time,sat,nav,dtss,&varc)) return 0;
    
    time_tt=timeadd(time,tt);
    if (!od_rtk_preceph_pephpos(time_tt,sat,nav,rst,dtst,NULL,NULL)||
        !od_rtk_preceph_pephclk(time_tt,sat,nav,dtst,NULL)) return 0;
    
    /* satellite antenna offset correction */
    if (opt) {
        satantoff(time,rss,sat,nav,dant);
    }
    for (i=0;i<3;i++) {
        rs[i  ]=rss[i]+dant[i];
        rs[i+3]=(rst[i]-rss[i])/tt;
    }
    /* relativistic effect correction */
    if (dtss[0]!=0.0) {
        dts[0]=dtss[0]-2.0*dot(rs,rs+3,3)/CLIGHT/CLIGHT;
        dts[1]=(dtst[0]-dtss[0])/tt;
    }
    else { /* no precise clock */
        dts[0]=dts[1]=0.0;
    }
    if (var) *var=vare+varc;
    
    return 1;
}


#pragma pop_macro("SQR")
#pragma pop_macro("NMAX")
#pragma pop_macro("MAXDTE")
#pragma pop_macro("EXTERR_EPH")
#pragma pop_macro("EXTERR_CLK")


/* ===== Embedded sbas.c ===== */
#pragma push_macro("WEEKOFFSET")
#undef WEEKOFFSET

/*------------------------------------------------------------------------------
* sbas.c : sbas functions
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* option : -DRRCENA  enable rrc correction
*          
* references :
*     [1] RTCA/DO-229C, Minimum operational performanc standards for global
*         positioning system/wide area augmentation system airborne equipment,
*         RTCA inc, November 28, 2001
*     [2] IS-QZSS v.1.1, Quasi-Zenith Satellite System Navigation Service
*         Interface Specification for QZSS, Japan Aerospace Exploration Agency,
*         July 31, 2009
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/10/14 1.0  new
*           2009/01/24 1.1  modify sbspntpos() api
*                           improve fast/ion correction update
*           2009/04/08 1.2  move function crc24q() to rcvlog.c
*                           support glonass, galileo and qzss
*           2009/06/08 1.3  modify sbsupdatestat()
*                           delete sbssatpos()
*           2009/12/12 1.4  support glonass
*           2010/01/22 1.5  support ems (egnos message service) format
*           2010/06/10 1.6  added api:
*                               sbssatcorr(),sbstropcorr(),sbsioncorr(),
*                               sbsupdatecorr()
*                           changed api:
*                               sbsreadmsgt(),sbsreadmsg()
*                           deleted api:
*                               sbspntpos(),sbsupdatestat()
*           2010/08/16 1.7  not reject udre==14 or give==15 correction message
*                           (2.4.0_p4)
*           2011/01/15 1.8  use api ionppp()
*                           add prn mask of qzss for qzss L1SAIF
*           2016/07/29 1.9  crc24q() -> rtk_crc24q()
*           2020/11/30 1.10 use integer types in stdint.h
*-----------------------------------------------------------------------------*/

/* constants -----------------------------------------------------------------*/

#define WEEKOFFSET  1024        /* gps week offset for NovAtel OEM-3 */

/* sbas igp definition -------------------------------------------------------*/
static const int16_t
x1[]={-75,-65,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,- 5,  0,  5, 10, 15, 20,
       25, 30, 35, 40, 45, 50, 55, 65, 75, 85},
x2[]={-55,-50,-45,-40,-35,-30,-25,-20,-15,-10, -5,  0,  5, 10, 15, 20, 25, 30,
       35, 40, 45, 50, 55},
x3[]={-75,-65,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,- 5,  0,  5, 10, 15, 20,
       25, 30, 35, 40, 45, 50, 55, 65, 75},
x4[]={-85,-75,-65,-55,-50,-45,-40,-35,-30,-25,-20,-15,-10,- 5,  0,  5, 10, 15,
       20, 25, 30, 35, 40, 45, 50, 55, 65, 75},
x5[]={-180,-175,-170,-165,-160,-155,-150,-145,-140,-135,-130,-125,-120,-115,
      -110,-105,-100,- 95,- 90,- 85,- 80,- 75,- 70,- 65,- 60,- 55,- 50,- 45,
      - 40,- 35,- 30,- 25,- 20,- 15,- 10,-  5,   0,   5,  10,  15,  20,  25,
        30,  35,  40,  45,  50,  55,  60,  65,  70,  75,  80,  85,  90,  95,
       100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165,
       170, 175},
x6[]={-180,-170,-160,-150,-140,-130,-120,-110,-100,- 90,- 80,- 70,- 60,- 50,
      - 40,- 30,- 20,- 10,   0,  10,  20,  30,  40,  50,  60,  70,  80,  90,
       100, 110, 120, 130, 140, 150, 160, 170},
x7[]={-180,-150,-120,- 90,- 60,- 30,   0,  30,  60,  90, 120, 150},
x8[]={-170,-140,-110,- 80,- 50,- 20,  10,  40,  70, 100, 130, 160};

EXPORT const sbsigpband_t igpband1[9][8]={ /* band 0-8 */
    {{-180,x1,  1, 28},{-175,x2, 29, 51},{-170,x3, 52, 78},{-165,x2, 79,101},
     {-160,x3,102,128},{-155,x2,129,151},{-150,x3,152,178},{-145,x2,179,201}},
    {{-140,x4,  1, 28},{-135,x2, 29, 51},{-130,x3, 52, 78},{-125,x2, 79,101},
     {-120,x3,102,128},{-115,x2,129,151},{-110,x3,152,178},{-105,x2,179,201}},
    {{-100,x3,  1, 27},{- 95,x2, 28, 50},{- 90,x1, 51, 78},{- 85,x2, 79,101},
     {- 80,x3,102,128},{- 75,x2,129,151},{- 70,x3,152,178},{- 65,x2,179,201}},
    {{- 60,x3,  1, 27},{- 55,x2, 28, 50},{- 50,x4, 51, 78},{- 45,x2, 79,101},
     {- 40,x3,102,128},{- 35,x2,129,151},{- 30,x3,152,178},{- 25,x2,179,201}},
    {{- 20,x3,  1, 27},{- 15,x2, 28, 50},{- 10,x3, 51, 77},{-  5,x2, 78,100},
     {   0,x1,101,128},{   5,x2,129,151},{  10,x3,152,178},{  15,x2,179,201}},
    {{  20,x3,  1, 27},{  25,x2, 28, 50},{  30,x3, 51, 77},{  35,x2, 78,100},
     {  40,x4,101,128},{  45,x2,129,151},{  50,x3,152,178},{  55,x2,179,201}},
    {{  60,x3,  1, 27},{  65,x2, 28, 50},{  70,x3, 51, 77},{  75,x2, 78,100},
     {  80,x3,101,127},{  85,x2,128,150},{  90,x1,151,178},{  95,x2,179,201}},
    {{ 100,x3,  1, 27},{ 105,x2, 28, 50},{ 110,x3, 51, 77},{ 115,x2, 78,100},
     { 120,x3,101,127},{ 125,x2,128,150},{ 130,x4,151,178},{ 135,x2,179,201}},
    {{ 140,x3,  1, 27},{ 145,x2, 28, 50},{ 150,x3, 51, 77},{ 155,x2, 78,100},
     { 160,x3,101,127},{ 165,x2,128,150},{ 170,x3,151,177},{ 175,x2,178,200}}
};
EXPORT const sbsigpband_t igpband2[2][5]={ /* band 9-10 */
    {{  60,x5,  1, 72},{  65,x6, 73,108},{  70,x6,109,144},{  75,x6,145,180},
     {  85,x7,181,192}},
    {{- 60,x5,  1, 72},{- 65,x6, 73,108},{- 70,x6,109,144},{- 75,x6,145,180},
     {- 85,x8,181,192}}
};
/* extract field from line ---------------------------------------------------*/
static char *od_rtk_sbas_getfield(char *p, int pos)
{
    for (pos--;pos>0;pos--,p++) if (!(p=strchr(p,','))) return NULL;
    return p;
}
/* variance of fast correction (udre=UDRE+1) ---------------------------------*/
static double od_rtk_sbas_varfcorr(int udre)
{
    const double var[14]={
        0.052,0.0924,0.1444,0.283,0.4678,0.8315,1.2992,1.8709,2.5465,3.326,
        5.1968,20.7870,230.9661,2078.695
    };
    return 0<udre&&udre<=14?var[udre-1]:0.0;
}
/* variance of ionosphere correction (give=GIVEI+1) --------------------------*/
static double od_rtk_sbas_varicorr(int give)
{
    const double var[15]={
        0.0084,0.0333,0.0749,0.1331,0.2079,0.2994,0.4075,0.5322,0.6735,0.8315,
        1.1974,1.8709,3.326,20.787,187.0826
    };
    return 0<give&&give<=15?var[give-1]:0.0;
}
/* fast correction degradation -----------------------------------------------*/
static double od_rtk_sbas_degfcorr(int ai)
{
    const double degf[16]={
        0.00000,0.00005,0.00009,0.00012,0.00015,0.00020,0.00030,0.00045,
        0.00060,0.00090,0.00150,0.00210,0.00270,0.00330,0.00460,0.00580
    };
    return 0<ai&&ai<=15?degf[ai]:0.0058;
}
/* decode type 1: prn masks --------------------------------------------------*/
static int od_rtk_sbas_decode_sbstype1(const sbsmsg_t *msg, sbssat_t *sbssat)
{
    int i,n,sat;
    
    trace(4,"decode_sbstype1:\n");
    
    for (i=1,n=0;i<=210&&n<MAXSAT;i++) {
        if (getbitu(msg->msg,13+i,1)) {
           if      (i<= 37) sat=satno(SYS_GPS,i);    /*   0- 37: gps */
           else if (i<= 61) sat=satno(SYS_GLO,i-37); /*  38- 61: glonass */
           else if (i<=119) sat=0;                   /*  62-119: future gnss */
           else if (i<=138) sat=satno(SYS_SBS,i);    /* 120-138: geo/waas */
           else if (i<=182) sat=0;                   /* 139-182: reserved */
           else if (i<=192) sat=satno(SYS_SBS,i+10); /* 183-192: qzss ref [2] */
           else if (i<=202) sat=satno(SYS_QZS,i);    /* 193-202: qzss ref [2] */
           else             sat=0;                   /* 203-   : reserved */
           sbssat->sat[n++].sat=sat;
        }
    }
    sbssat->iodp=getbitu(msg->msg,224,2);
    sbssat->nsat=n;
    
    trace(5,"decode_sbstype1: nprn=%d iodp=%d\n",n,sbssat->iodp);
    return 1;
}
/* decode type 2-5,0: fast corrections ---------------------------------------*/
static int od_rtk_sbas_decode_sbstype2(const sbsmsg_t *msg, sbssat_t *sbssat)
{
    int i,j,iodf,type,udre;
    double prc,dt;
    gtime_t t0;
    
    trace(4,"decode_sbstype2:\n");
    
    if (sbssat->iodp!=(int)getbitu(msg->msg,16,2)) return 0;
    
    type=getbitu(msg->msg, 8,6);
    iodf=getbitu(msg->msg,14,2);
    
    for (i=0;i<13;i++) {
        if ((j=13*((type==0?2:type)-2)+i)>=sbssat->nsat) break;
        udre=getbitu(msg->msg,174+4*i,4);
        t0 =sbssat->sat[j].fcorr.t0;
        prc=sbssat->sat[j].fcorr.prc;
        sbssat->sat[j].fcorr.t0=gpst2time(msg->week,msg->tow);
        sbssat->sat[j].fcorr.prc=getbits(msg->msg,18+i*12,12)*0.125f;
        sbssat->sat[j].fcorr.udre=udre+1;
        dt=timediff(sbssat->sat[j].fcorr.t0,t0);
        if (t0.time==0||dt<=0.0||18.0<dt||sbssat->sat[j].fcorr.ai==0) {
            sbssat->sat[j].fcorr.rrc=0.0;
            sbssat->sat[j].fcorr.dt=0.0;
        }
        else {
            sbssat->sat[j].fcorr.rrc=(sbssat->sat[j].fcorr.prc-prc)/dt;
            sbssat->sat[j].fcorr.dt=dt;
        }
        sbssat->sat[j].fcorr.iodf=iodf;
    }
    trace(5,"decode_sbstype2: type=%d iodf=%d\n",type,iodf);
    return 1;
}
/* decode type 6: integrity info ---------------------------------------------*/
static int od_rtk_sbas_decode_sbstype6(const sbsmsg_t *msg, sbssat_t *sbssat)
{
    int i,iodf[4],udre;
    
    trace(4,"decode_sbstype6:\n");
    
    for (i=0;i<4;i++) {
        iodf[i]=getbitu(msg->msg,14+i*2,2);
    }
    for (i=0;i<sbssat->nsat&&i<MAXSAT;i++) {
        if (sbssat->sat[i].fcorr.iodf!=iodf[i/13]) continue;
        udre=getbitu(msg->msg,22+i*4,4);
        sbssat->sat[i].fcorr.udre=udre+1;
    }
    trace(5,"decode_sbstype6: iodf=%d %d %d %d\n",iodf[0],iodf[1],iodf[2],iodf[3]);
    return 1;
}
/* decode type 7: fast correction degradation factor -------------------------*/
static int od_rtk_sbas_decode_sbstype7(const sbsmsg_t *msg, sbssat_t *sbssat)
{
    int i;
    
    trace(4,"decode_sbstype7\n");
    
    if (sbssat->iodp!=(int)getbitu(msg->msg,18,2)) return 0;
    
    sbssat->tlat=getbitu(msg->msg,14,4);
    
    for (i=0;i<sbssat->nsat&&i<MAXSAT;i++) {
        sbssat->sat[i].fcorr.ai=getbitu(msg->msg,22+i*4,4);
    }
    return 1;
}
/* decode type 9: geo navigation message -------------------------------------*/
static int od_rtk_sbas_decode_sbstype9(const sbsmsg_t *msg, nav_t *nav)
{
    seph_t seph={0};
    int i,sat,t;
    
    trace(4,"decode_sbstype9:\n");
    
    if (!(sat=satno(SYS_SBS,msg->prn))) {
        trace(2,"invalid prn in sbas type 9: prn=%3d\n",msg->prn);
        return 0;
    }
    t=(int)getbitu(msg->msg,22,13)*16-(int)msg->tow%86400;
    if      (t<=-43200) t+=86400;
    else if (t>  43200) t-=86400;
    seph.sat=sat;
    seph.t0 =gpst2time(msg->week,msg->tow+t);
    seph.tof=gpst2time(msg->week,msg->tow);
    seph.sva=getbitu(msg->msg,35,4);
    seph.svh=seph.sva==15?1:0; /* unhealthy if ura==15 */
    
    seph.pos[0]=getbits(msg->msg, 39,30)*0.08;
    seph.pos[1]=getbits(msg->msg, 69,30)*0.08;
    seph.pos[2]=getbits(msg->msg, 99,25)*0.4;
    seph.vel[0]=getbits(msg->msg,124,17)*0.000625;
    seph.vel[1]=getbits(msg->msg,141,17)*0.000625;
    seph.vel[2]=getbits(msg->msg,158,18)*0.004;
    seph.acc[0]=getbits(msg->msg,176,10)*0.0000125;
    seph.acc[1]=getbits(msg->msg,186,10)*0.0000125;
    seph.acc[2]=getbits(msg->msg,196,10)*0.0000625;
    
    seph.af0=getbits(msg->msg,206,12)*P2_31;
    seph.af1=getbits(msg->msg,218, 8)*P2_39/2.0;
    
    i=msg->prn-MINPRNSBS;
    if (!nav->seph||fabs(timediff(nav->seph[i].t0,seph.t0))<1E-3) { /* not change */
        return 0;
    }
    nav->seph[NSATSBS+i]=nav->seph[i]; /* previous */
    nav->seph[i]=seph;                 /* current */

    trace(5,"decode_sbstype9: prn=%d\n",msg->prn);
    return 1;
}
/* decode type 18: ionospheric grid point masks ------------------------------*/
static int od_rtk_sbas_decode_sbstype18(const sbsmsg_t *msg, sbsion_t *sbsion)
{
    const sbsigpband_t *p;
    int i,j,n,m,band=getbitu(msg->msg,18,4);
    
    trace(4,"decode_sbstype18:\n");
    
    if      (0<=band&&band<= 8) {p=igpband1[band  ]; m=8;}
    else if (9<=band&&band<=10) {p=igpband2[band-9]; m=5;}
    else return 0;
    
    sbsion[band].iodi=(int16_t)getbitu(msg->msg,22,2);
    
    for (i=1,n=0;i<=201;i++) {
        if (!getbitu(msg->msg,23+i,1)) continue;
        for (j=0;j<m;j++) {
            if (i<p[j].bits||p[j].bite<i) continue;
            sbsion[band].igp[n].lat=band<=8?p[j].y[i-p[j].bits]:p[j].x;
            sbsion[band].igp[n++].lon=band<=8?p[j].x:p[j].y[i-p[j].bits];
            break;
        }
    }
    sbsion[band].nigp=n;
    
    trace(5,"decode_sbstype18: band=%d nigp=%d\n",band,n);
    return 1;
}
/* decode half long term correction (vel code=0) -----------------------------*/
static int od_rtk_sbas_decode_longcorr0(const sbsmsg_t *msg, int p, sbssat_t *sbssat)
{
    int i,n=getbitu(msg->msg,p,6);
    
    trace(4,"decode_longcorr0:\n");
    
    if (n==0||n>MAXSAT) return 0;
    
    sbssat->sat[n-1].lcorr.iode=getbitu(msg->msg,p+6,8);
    
    for (i=0;i<3;i++) {
        sbssat->sat[n-1].lcorr.dpos[i]=getbits(msg->msg,p+14+9*i,9)*0.125;
        sbssat->sat[n-1].lcorr.dvel[i]=0.0;
    }
    sbssat->sat[n-1].lcorr.daf0=getbits(msg->msg,p+41,10)*P2_31;
    sbssat->sat[n-1].lcorr.daf1=0.0;
    sbssat->sat[n-1].lcorr.t0=gpst2time(msg->week,msg->tow);
    
    trace(5,"decode_longcorr0:sat=%2d\n",sbssat->sat[n-1].sat);
    return 1;
}
/* decode half long term correction (vel code=1) -----------------------------*/
static int od_rtk_sbas_decode_longcorr1(const sbsmsg_t *msg, int p, sbssat_t *sbssat)
{
    int i,n=getbitu(msg->msg,p,6),t;
    
    trace(4,"decode_longcorr1:\n");
    
    if (n==0||n>MAXSAT) return 0;
    
    sbssat->sat[n-1].lcorr.iode=getbitu(msg->msg,p+6,8);
    
    for (i=0;i<3;i++) {
        sbssat->sat[n-1].lcorr.dpos[i]=getbits(msg->msg,p+14+i*11,11)*0.125;
        sbssat->sat[n-1].lcorr.dvel[i]=getbits(msg->msg,p+58+i* 8, 8)*P2_11;
    }
    sbssat->sat[n-1].lcorr.daf0=getbits(msg->msg,p+47,11)*P2_31;
    sbssat->sat[n-1].lcorr.daf1=getbits(msg->msg,p+82, 8)*P2_39;
    t=(int)getbitu(msg->msg,p+90,13)*16-(int)msg->tow%86400;
    if      (t<=-43200) t+=86400;
    else if (t>  43200) t-=86400;
    sbssat->sat[n-1].lcorr.t0=gpst2time(msg->week,msg->tow+t);
    
    trace(5,"decode_longcorr1: sat=%2d\n",sbssat->sat[n-1].sat);
    return 1;
}
/* decode half long term correction ------------------------------------------*/
static int od_rtk_sbas_decode_longcorrh(const sbsmsg_t *msg, int p, sbssat_t *sbssat)
{
    trace(4,"decode_longcorrh:\n");
    
    if (getbitu(msg->msg,p,1)==0) { /* vel code=0 */
        if (sbssat->iodp==(int)getbitu(msg->msg,p+103,2)) {
            return od_rtk_sbas_decode_longcorr0(msg,p+ 1,sbssat)&&
                   od_rtk_sbas_decode_longcorr0(msg,p+52,sbssat);
        }
    }
    else if (sbssat->iodp==(int)getbitu(msg->msg,p+104,2)) {
        return od_rtk_sbas_decode_longcorr1(msg,p+1,sbssat);
    }
    return 0;
}
/* decode type 24: mixed fast/long term correction ---------------------------*/
static int od_rtk_sbas_decode_sbstype24(const sbsmsg_t *msg, sbssat_t *sbssat)
{
    int i,j,iodf,blk,udre;
    
    trace(4,"decode_sbstype24:\n");
    
    if (sbssat->iodp!=(int)getbitu(msg->msg,110,2)) return 0; /* check IODP */
    
    blk =getbitu(msg->msg,112,2);
    iodf=getbitu(msg->msg,114,2);
    
    for (i=0;i<6;i++) {
        if ((j=13*blk+i)>=sbssat->nsat) break;
        udre=getbitu(msg->msg,86+4*i,4);
        
        sbssat->sat[j].fcorr.t0  =gpst2time(msg->week,msg->tow);
        sbssat->sat[j].fcorr.prc =getbits(msg->msg,14+i*12,12)*0.125f;
        sbssat->sat[j].fcorr.udre=udre+1;
        sbssat->sat[j].fcorr.iodf=iodf;
    }
    return od_rtk_sbas_decode_longcorrh(msg,120,sbssat);
}
/* decode type 25: long term satellite error correction ----------------------*/
static int od_rtk_sbas_decode_sbstype25(const sbsmsg_t *msg, sbssat_t *sbssat)
{
    trace(4,"decode_sbstype25:\n");
    
    return od_rtk_sbas_decode_longcorrh(msg,14,sbssat)&&od_rtk_sbas_decode_longcorrh(msg,120,sbssat);
}
/* decode type 26: ionospheric deley corrections -----------------------------*/
static int od_rtk_sbas_decode_sbstype26(const sbsmsg_t *msg, sbsion_t *sbsion)
{
    int i,j,block,delay,give,band=getbitu(msg->msg,14,4);
    
    trace(4,"decode_sbstype26:\n");
    
    if (band>MAXBAND||sbsion[band].iodi!=(int)getbitu(msg->msg,217,2)) return 0;
    
    block=getbitu(msg->msg,18,4);
    
    for (i=0;i<15;i++) {
        if ((j=block*15+i)>=sbsion[band].nigp) continue;
        give=getbitu(msg->msg,22+i*13+9,4);
        
        delay=getbitu(msg->msg,22+i*13,9);
        sbsion[band].igp[j].t0=gpst2time(msg->week,msg->tow);
        sbsion[band].igp[j].delay=delay==0x1FF?0.0f:delay*0.125f;
        sbsion[band].igp[j].give=give+1;
        
        if (sbsion[band].igp[j].give>=16) {
            sbsion[band].igp[j].give=0;
        }
    }
    trace(5,"decode_sbstype26: band=%d block=%d\n",band,block);
    return 1;
}
/* update sbas corrections -----------------------------------------------------
* update sbas correction parameters in navigation data with a sbas message
* args   : sbsmg_t  *msg    I   sbas message
*          nav_t    *nav    IO  navigation data
* return : message type (-1: error or not supported type)
* notes  : nav->seph must point to seph[NSATSBS*2] (array of seph_t)
*               seph[prn-MINPRNSBS+1]          : sat prn current epehmeris 
*               seph[prn-MINPRNSBS+1+MAXPRNSBS]: sat prn previous epehmeris 
*-----------------------------------------------------------------------------*/
extern int sbsupdatecorr(const sbsmsg_t *msg, nav_t *nav)
{
    int type=getbitu(msg->msg,8,6),stat=-1;
    
    trace(3,"sbsupdatecorr: type=%d\n",type);
    
    if (msg->week==0) return -1;
    
    switch (type) {
        case  0: stat=od_rtk_sbas_decode_sbstype2 (msg,&nav->sbssat); break;
        case  1: stat=od_rtk_sbas_decode_sbstype1 (msg,&nav->sbssat); break;
        case  2:
        case  3:
        case  4:
        case  5: stat=od_rtk_sbas_decode_sbstype2 (msg,&nav->sbssat); break;
        case  6: stat=od_rtk_sbas_decode_sbstype6 (msg,&nav->sbssat); break;
        case  7: stat=od_rtk_sbas_decode_sbstype7 (msg,&nav->sbssat); break;
        case  9: stat=od_rtk_sbas_decode_sbstype9 (msg,nav);          break;
        case 18: stat=od_rtk_sbas_decode_sbstype18(msg,nav ->sbsion); break;
        case 24: stat=od_rtk_sbas_decode_sbstype24(msg,&nav->sbssat); break;
        case 25: stat=od_rtk_sbas_decode_sbstype25(msg,&nav->sbssat); break;
        case 26: stat=od_rtk_sbas_decode_sbstype26(msg,nav ->sbsion); break;
        case 63: break; /* null message */
        
        /*default: trace(2,"unsupported sbas message: type=%d\n",type); break;*/
    }
    return stat?type:-1;
}
/* read sbas log file --------------------------------------------------------*/
static void od_rtk_sbas_readmsgs(const char *file, int sel, gtime_t ts, gtime_t te,
                     sbs_t *sbs)
{
    sbsmsg_t *sbs_msgs;
    int i,week,prn,ch,msg;
    uint32_t b;
    double tow,ep[6]={0};
    char buff[256],*p;
    gtime_t time;
    FILE *fp;
    
    trace(3,"readmsgs: file=%s sel=%d\n",file,sel);
    
    if (!(fp=fopen(file,"r"))) {
        trace(2,"sbas message file open error: %s\n",file);
        return;
    }
    while (fgets(buff,sizeof(buff),fp)) {
        if (sscanf(buff,"%d %lf %d",&week,&tow,&prn)==3&&(p=strstr(buff,": "))) {
            p+=2; /* rtklib form */
        }
        else if (sscanf(buff,"%d %lf %lf %lf %lf %lf %lf %d",
                        &prn,ep,ep+1,ep+2,ep+3,ep+4,ep+5,&msg)==8) {
            /* ems (EGNOS Message Service) form */
            ep[0]+=ep[0]<70.0?2000.0:1900.0;
            tow=time2gpst(epoch2time(ep),&week);
            p=buff+(msg>=10?25:24);
        }
        else if (!strncmp(buff,"#RAWWAASFRAMEA",14)) { /* NovAtel OEM4/V */
            if (!(p=od_rtk_sbas_getfield(buff,6))) continue;
            if (sscanf(p,"%d,%lf",&week,&tow)<2) continue;
            if (!(p=strchr(++p,';'))) continue;
            if (sscanf(++p,"%d,%d",&ch,&prn)<2) continue;
            if (!(p=od_rtk_sbas_getfield(p,4))) continue;
        }
        else if (!strncmp(buff,"$FRMA",5)) { /* NovAtel OEM3 */
            if (!(p=od_rtk_sbas_getfield(buff,2))) continue;
            if (sscanf(p,"%d,%lf,%d",&week,&tow,&prn)<3) continue;
            if (!(p=od_rtk_sbas_getfield(p,6))) continue;
            if (week<WEEKOFFSET) week+=WEEKOFFSET;
        }
        else continue;
        
        if (sel!=0&&sel!=prn) continue;
        
        time=gpst2time(week,tow);
        
        if (!screent(time,ts,te,0.0)) continue;
        
        if (sbs->n>=sbs->nmax) {
            sbs->nmax=sbs->nmax==0?1024:sbs->nmax*2;
            if (!(sbs_msgs=(sbsmsg_t *)realloc(sbs->msgs,sbs->nmax*sizeof(sbsmsg_t)))) {
                trace(1,"readsbsmsg malloc error: nmax=%d\n",sbs->nmax);
                free(sbs->msgs); sbs->msgs=NULL; sbs->n=sbs->nmax=0;
                return;
            }
            sbs->msgs=sbs_msgs;
        }
        sbs->msgs[sbs->n].week=week;
        sbs->msgs[sbs->n].tow=(int)(tow+0.5);
        sbs->msgs[sbs->n].prn=prn;
        for (i=0;i<29;i++) sbs->msgs[sbs->n].msg[i]=0;
        for (i=0;*(p-1)&&*p&&i<29;p+=2,i++) {
            if (sscanf(p,"%2X",&b)==1) sbs->msgs[sbs->n].msg[i]=(uint8_t)b;
        }
        sbs->msgs[sbs->n++].msg[28]&=0xC0;
    }
    fclose(fp);
}
/* compare sbas messages -----------------------------------------------------*/
static int od_rtk_sbas_cmpmsgs(const void *p1, const void *p2)
{
    sbsmsg_t *q1=(sbsmsg_t *)p1,*q2=(sbsmsg_t *)p2;
    return q1->week!=q2->week?q1->week-q2->week:
           (q1->tow<q2->tow?-1:(q1->tow>q2->tow?1:q1->prn-q2->prn));
}
/* read sbas message file ------------------------------------------------------
* read sbas message file
* args   : char     *file   I   sbas message file (wind-card * is expanded)
*          int      sel     I   sbas satellite prn number selection (0:all)
*         (gtime_t  ts      I   start time)
*         (gtime_t  te      I   end time  )
*          sbs_t    *sbs    IO  sbas messages
* return : number of sbas messages
* notes  : sbas message are appended and sorted. before calling the funciton, 
*          sbs->n, sbs->nmax and sbs->msgs must be set properly. (initially
*          sbs->n=sbs->nmax=0, sbs->msgs=NULL)
*          only the following file extentions after wild card expanded are valid
*          to read. others are skipped
*          .sbs, .SBS, .ems, .EMS
*-----------------------------------------------------------------------------*/
extern int sbsreadmsgt(const char *file, int sel, gtime_t ts, gtime_t te,
                       sbs_t *sbs)
{
    char *efiles[MAXEXFILE]={0},*ext;
    int i,n;
    
    trace(3,"sbsreadmsgt: file=%s sel=%d\n",file,sel);
    
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return 0;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
        if (!(ext=strrchr(efiles[i],'.'))) continue;
        if (strcmp(ext,".sbs")&&strcmp(ext,".SBS")&&
            strcmp(ext,".ems")&&strcmp(ext,".EMS")) continue;
        
        od_rtk_sbas_readmsgs(efiles[i],sel,ts,te,sbs);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    /* sort messages */
    if (sbs->n>0) {
        qsort(sbs->msgs,sbs->n,sizeof(sbsmsg_t),od_rtk_sbas_cmpmsgs);
    }
    return sbs->n;
}
extern int sbsreadmsg(const char *file, int sel, sbs_t *sbs)
{
    gtime_t ts={0},te={0};
    
    trace(3,"sbsreadmsg: file=%s sel=%d\n",file,sel);
    
    return sbsreadmsgt(file,sel,ts,te,sbs);
}
/* output sbas messages --------------------------------------------------------
* output sbas message record to output file in rtklib sbas log format
* args   : FILE   *fp       I   output file pointer
*          sbsmsg_t *sbsmsg I   sbas messages
* return : none
*-----------------------------------------------------------------------------*/
extern void sbsoutmsg(FILE *fp, sbsmsg_t *sbsmsg)
{
    int i,prn=sbsmsg->prn,type=sbsmsg->msg[1]>>2;
    
    trace(4,"sbsoutmsg:\n");
    
    fprintf(fp,"%4d %6d %3d %2d : ",sbsmsg->week,sbsmsg->tow,prn,type);
    for (i=0;i<29;i++) fprintf(fp,"%02X",sbsmsg->msg[i]);
    fprintf(fp,"\n");
}
/* search igps ---------------------------------------------------------------*/
static void od_rtk_sbas_searchigp(gtime_t time, const double *pos, const sbsion_t *ion,
                      const sbsigp_t **igp, double *x, double *y)
{
    int i,latp[2],lonp[4];
    double lat=pos[0]*R2D,lon=pos[1]*R2D;
    const sbsigp_t *p;
    
    trace(4,"searchigp: pos=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D);
    
    if (lon>=180.0) lon-=360.0;
    if (-55.0<=lat&&lat<55.0) {
        latp[0]=(int)floor(lat/5.0)*5;
        latp[1]=latp[0]+5;
        lonp[0]=lonp[1]=(int)floor(lon/5.0)*5;
        lonp[2]=lonp[3]=lonp[0]+5;
        *x=(lon-lonp[0])/5.0;
        *y=(lat-latp[0])/5.0;
    }
    else {
        latp[0]=(int)floor((lat-5.0)/10.0)*10+5;
        latp[1]=latp[0]+10;
        lonp[0]=lonp[1]=(int)floor(lon/10.0)*10;
        lonp[2]=lonp[3]=lonp[0]+10;
        *x=(lon-lonp[0])/10.0;
        *y=(lat-latp[0])/10.0;
        if (75.0<=lat&&lat<85.0) {
            lonp[1]=(int)floor(lon/90.0)*90;
            lonp[3]=lonp[1]+90;
        }
        else if (-85.0<=lat&&lat<-75.0) {
            lonp[0]=(int)floor((lon-50.0)/90.0)*90+40;
            lonp[2]=lonp[0]+90;
        }
        else if (lat>=85.0) {
            for (i=0;i<4;i++) lonp[i]=(int)floor(lon/90.0)*90;
        }
        else if (lat<-85.0) {
            for (i=0;i<4;i++) lonp[i]=(int)floor((lon-50.0)/90.0)*90+40;
        }
    }
    for (i=0;i<4;i++) if (lonp[i]==180) lonp[i]=-180;
    for (i=0;i<=MAXBAND;i++) {
        for (p=ion[i].igp;p<ion[i].igp+ion[i].nigp;p++) {
            if (p->t0.time==0) continue;
            if      (p->lat==latp[0]&&p->lon==lonp[0]&&p->give>0) igp[0]=p;
            else if (p->lat==latp[1]&&p->lon==lonp[1]&&p->give>0) igp[1]=p;
            else if (p->lat==latp[0]&&p->lon==lonp[2]&&p->give>0) igp[2]=p;
            else if (p->lat==latp[1]&&p->lon==lonp[3]&&p->give>0) igp[3]=p;
            if (igp[0]&&igp[1]&&igp[2]&&igp[3]) return;
        }
    }
}
/* sbas ionospheric delay correction -------------------------------------------
* compute sbas ionosphric delay correction
* args   : gtime_t  time    I   time
*          nav_t    *nav    I   navigation data
*          double   *pos    I   receiver position {lat,lon,height} (rad/m)
*          double   *azel   I   satellite azimuth/elavation angle (rad)
*          double   *delay  O   slant ionospheric delay (L1) (m)
*          double   *var    O   variance of ionospheric delay (m^2)
* return : status (1:ok, 0:no correction)
* notes  : before calling the function, sbas ionosphere correction parameters
*          in navigation data (nav->sbsion) must be set by callig 
*          sbsupdatecorr()
*-----------------------------------------------------------------------------*/
extern int sbsioncorr(gtime_t time, const nav_t *nav, const double *pos,
                      const double *azel, double *delay, double *var)
{
    const double re=6378.1363,hion=350.0;
    int i,err=0;
    double fp,posp[2],x=0.0,y=0.0,t,w[4]={0};
    const sbsigp_t *igp[4]={0}; /* {ws,wn,es,en} */
    
    trace(4,"sbsioncorr: pos=%.3f %.3f azel=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D,
          azel[0]*R2D,azel[1]*R2D);
    
    *delay=*var=0.0;
    if (pos[2]<-100.0||azel[1]<=0) return 1;
    
    /* ipp (ionospheric pierce point) position */
    fp=ionppp(pos,azel,re,hion,posp);
    
    /* search igps around ipp */
    od_rtk_sbas_searchigp(time,posp,nav->sbsion,igp,&x,&y);
    
    /* weight of igps */
    if (igp[0]&&igp[1]&&igp[2]&&igp[3]) {
        w[0]=(1.0-x)*(1.0-y); w[1]=(1.0-x)*y; w[2]=x*(1.0-y); w[3]=x*y;
    }
    else if (igp[0]&&igp[1]&&igp[2]) {
        w[1]=y; w[2]=x;
        if ((w[0]=1.0-w[1]-w[2])<0.0) err=1;
    }
    else if (igp[0]&&igp[2]&&igp[3]) {
        w[0]=1.0-x; w[3]=y;
        if ((w[2]=1.0-w[0]-w[3])<0.0) err=1;
    }
    else if (igp[0]&&igp[1]&&igp[3]) {
        w[0]=1.0-y; w[3]=x;
        if ((w[1]=1.0-w[0]-w[3])<0.0) err=1;
    }
    else if (igp[1]&&igp[2]&&igp[3]) {
        w[1]=1.0-x; w[2]=1.0-y;
        if ((w[3]=1.0-w[1]-w[2])<0.0) err=1;
    }
    else err=1;
    
    if (err) {
        trace(2,"no sbas iono correction: lat=%3.0f lon=%4.0f\n",posp[0]*R2D,
              posp[1]*R2D);
        return 0;
    }
    for (i=0;i<4;i++) {
        if (!igp[i]) continue;
        t=timediff(time,igp[i]->t0);
        *delay+=w[i]*igp[i]->delay;
        *var+=w[i]*od_rtk_sbas_varicorr(igp[i]->give)*9E-8*fabs(t);
    }
    *delay*=fp; *var*=fp*fp;
    
    trace(5,"sbsioncorr: dion=%7.2f sig=%7.2f\n",*delay,sqrt(*var));
    return 1;
}
/* get meterological parameters ----------------------------------------------*/
static void od_rtk_sbas_getmet(double lat, double *met)
{
    static const double metprm[][10]={ /* lat=15,30,45,60,75 */
        {1013.25,299.65,26.31,6.30E-3,2.77,  0.00, 0.00,0.00,0.00E-3,0.00},
        {1017.25,294.15,21.79,6.05E-3,3.15, -3.75, 7.00,8.85,0.25E-3,0.33},
        {1015.75,283.15,11.66,5.58E-3,2.57, -2.25,11.00,7.24,0.32E-3,0.46},
        {1011.75,272.15, 6.78,5.39E-3,1.81, -1.75,15.00,5.36,0.81E-3,0.74},
        {1013.00,263.65, 4.11,4.53E-3,1.55, -0.50,14.50,3.39,0.62E-3,0.30}
    };
    int i,j;
    double a;
    lat=fabs(lat);
    if      (lat<=15.0) for (i=0;i<10;i++) met[i]=metprm[0][i];
    else if (lat>=75.0) for (i=0;i<10;i++) met[i]=metprm[4][i];
    else {
        j=(int)(lat/15.0); a=(lat-j*15.0)/15.0;
        for (i=0;i<10;i++) met[i]=(1.0-a)*metprm[j-1][i]+a*metprm[j][i];
    }
}
/* tropospheric delay correction -----------------------------------------------
* compute sbas tropospheric delay correction (mops model)
* args   : gtime_t time     I   time
*          double   *pos    I   receiver position {lat,lon,height} (rad/m)
*          double   *azel   I   satellite azimuth/elavation (rad)
*          double   *var    O   variance of troposphric error (m^2)
* return : slant tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double sbstropcorr(gtime_t time, const double *pos, const double *azel,
                          double *var)
{
    const double k1=77.604,k2=382000.0,rd=287.054,gm=9.784,g=9.80665;
    static double pos_[3]={0},zh=0.0,zw=0.0;
    int i;
    double c,met[10],sinel=sin(azel[1]),h=pos[2],m;
    
    trace(4,"sbstropcorr: pos=%.3f %.3f azel=%.3f %.3f\n",pos[0]*R2D,pos[1]*R2D,
          azel[0]*R2D,azel[1]*R2D);
    
    if (pos[2]<-100.0||10000.0<pos[2]||azel[1]<=0) {
        *var=0.0;
        return 0.0;
    }
    if (zh==0.0||fabs(pos[0]-pos_[0])>1E-7||fabs(pos[1]-pos_[1])>1E-7||
        fabs(pos[2]-pos_[2])>1.0) {
        od_rtk_sbas_getmet(pos[0]*R2D,met);
        c=cos(2.0*PI*(time2doy(time)-(pos[0]>=0.0?28.0:211.0))/365.25);
        for (i=0;i<5;i++) met[i]-=met[i+5]*c;
        zh=1E-6*k1*rd*met[0]/gm;
        zw=1E-6*k2*rd/(gm*(met[4]+1.0)-met[3]*rd)*met[2]/met[1];
        zh*=pow(1.0-met[3]*h/met[1],g/(rd*met[3]));
        zw*=pow(1.0-met[3]*h/met[1],(met[4]+1.0)*g/(rd*met[3])-1.0);
        for (i=0;i<3;i++) pos_[i]=pos[i];
    }
    m=1.001/sqrt(0.002001+sinel*sinel);
    *var=0.12*0.12*m*m;
    return (zh+zw)*m;
}
/* long term correction ------------------------------------------------------*/
static int od_rtk_sbas_sbslongcorr(gtime_t time, int sat, const sbssat_t *sbssat,
                       double *drs, double *ddts)
{
    const sbssatp_t *p;
    double t;
    int i;
    
    trace(3,"sbslongcorr: sat=%2d\n",sat);
    
    for (p=sbssat->sat;p<sbssat->sat+sbssat->nsat;p++) {
        if (p->sat!=sat||p->lcorr.t0.time==0) continue;
        t=timediff(time,p->lcorr.t0);
        if (fabs(t)>MAXSBSAGEL) {
            trace(2,"sbas long-term correction expired: %s sat=%2d t=%5.0f\n",
                  time_str(time,0),sat,t);
            return 0;
        }
        for (i=0;i<3;i++) drs[i]=p->lcorr.dpos[i]+p->lcorr.dvel[i]*t;
        *ddts=p->lcorr.daf0+p->lcorr.daf1*t;
        
        trace(5,"sbslongcorr: sat=%2d drs=%7.2f%7.2f%7.2f ddts=%7.2f\n",
              sat,drs[0],drs[1],drs[2],*ddts*CLIGHT);
        
        return 1;
    }
    /* if sbas satellite without correction, no correction applied */
    if (satsys(sat,NULL)==SYS_SBS) return 1;
    
    trace(2,"no sbas long-term correction: %s sat=%2d\n",time_str(time,0),sat);
    return 0;
}
/* fast correction -----------------------------------------------------------*/
static int od_rtk_sbas_sbsfastcorr(gtime_t time, int sat, const sbssat_t *sbssat,
                       double *prc, double *var)
{
    const sbssatp_t *p;
    double t;
    
    trace(3,"sbsfastcorr: sat=%2d\n",sat);
    
    for (p=sbssat->sat;p<sbssat->sat+sbssat->nsat;p++) {
        if (p->sat!=sat) continue;
        if (p->fcorr.t0.time==0) break;
        t=timediff(time,p->fcorr.t0)+sbssat->tlat;
        
        /* expire age of correction or UDRE==14 (not monitored) */
        if (fabs(t)>MAXSBSAGEF||p->fcorr.udre>=15) continue;
        *prc=p->fcorr.prc;
#ifdef RRCENA
        if (p->fcorr.ai>0&&fabs(t)<=8.0*p->fcorr.dt) {
            *prc+=p->fcorr.rrc*t;
        }
#endif
        *var=od_rtk_sbas_varfcorr(p->fcorr.udre)+od_rtk_sbas_degfcorr(p->fcorr.ai)*t*t/2.0;
        
        trace(5,"sbsfastcorr: sat=%3d prc=%7.2f sig=%7.2f t=%5.0f\n",sat,
              *prc,sqrt(*var),t);
        return 1;
    }
    trace(2,"no sbas fast correction: %s sat=%2d\n",time_str(time,0),sat);
    return 0;
}
/* sbas satellite ephemeris and clock correction -------------------------------
* correct satellite position and clock bias with sbas satellite corrections
* args   : gtime_t time     I   reception time
*          int    sat       I   satellite
*          nav_t  *nav      I   navigation data
*          double *rs       IO  sat position and corrected {x,y,z} (ecef) (m)
*          double *dts      IO  sat clock bias and corrected (s)
*          double *var      O   sat position and clock variance (m^2)
* return : status (1:ok,0:no correction)
* notes  : before calling the function, sbas satellite correction parameters 
*          in navigation data (nav->sbssat) must be set by callig
*          sbsupdatecorr().
*          satellite clock correction include long-term correction and fast
*          correction.
*          sbas clock correction is usually based on L1C/A code. TGD or DCB has
*          to be considered for other codes
*-----------------------------------------------------------------------------*/
extern int sbssatcorr(gtime_t time, int sat, const nav_t *nav, double *rs,
                      double *dts, double *var)
{
    double drs[3]={0},dclk=0.0,prc=0.0;
    int i;
    
    trace(3,"sbssatcorr : sat=%2d\n",sat);
    
    /* sbas long term corrections */
    if (!od_rtk_sbas_sbslongcorr(time,sat,&nav->sbssat,drs,&dclk)) {
        return 0;
    }
    /* sbas fast corrections */
    if (!od_rtk_sbas_sbsfastcorr(time,sat,&nav->sbssat,&prc,var)) {
        return 0;
    }
    for (i=0;i<3;i++) rs[i]+=drs[i];
    
    dts[0]+=dclk+prc/CLIGHT;
    
    trace(5,"sbssatcorr: sat=%2d drs=%6.3f %6.3f %6.3f dclk=%.3f %.3f var=%.3f\n",
          sat,drs[0],drs[1],drs[2],dclk,prc/CLIGHT,*var);
    
    return 1;
}
/* decode sbas message ---------------------------------------------------------
* decode sbas message frame words and check crc
* args   : gtime_t time     I   reception time
*          int    prn       I   sbas satellite prn number
*          uint32_t *word   I   message frame words (24bit x 10)
*          sbsmsg_t *sbsmsg O   sbas message
* return : status (1:ok,0:crc error)
*-----------------------------------------------------------------------------*/
extern int sbsdecodemsg(gtime_t time, int prn, const uint32_t *words,
                        sbsmsg_t *sbsmsg)
{
    int i,j;
    uint8_t f[29];
    double tow;
    
    trace(5,"sbsdecodemsg: prn=%d\n",prn);
    
    if (time.time==0) return 0;
    tow=time2gpst(time,&sbsmsg->week);
    sbsmsg->tow=(int)(tow+DTTOL);
    sbsmsg->prn=prn;
    for (i=0;i<7;i++) for (j=0;j<4;j++) {
        sbsmsg->msg[i*4+j]=(uint8_t)(words[i]>>((3-j)*8));
    }
    sbsmsg->msg[28]=(uint8_t)(words[7]>>18)&0xC0;
    for (i=28;i>0;i--) f[i]=(sbsmsg->msg[i]>>6)+(sbsmsg->msg[i-1]<<2);
    f[0]=sbsmsg->msg[0]>>6;
    
    return rtk_crc24q(f,29)==(words[7]&0xFFFFFF); /* check crc */
}


#pragma pop_macro("WEEKOFFSET")

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
/* ===== Pure C input processing and public API ===== */
#define OD_DAY_MS INT64_C(86400000)
#define OD_MAX_TIME INT64_C(253402300799000)
#define OD_MAX_CHUNK (1024U * 1024U)
#define OD_MAX_SENTENCE 1024U

typedef struct { int year, month, day; } od_date_t;
typedef struct {
    rtcm_t decoder;
    prcopt_t options;
    int has_time;
} od_rtcm_solver_t;

struct od_context {
    od_config_t config;
    od_fit_observation_t *observations;
    od_fit_context_t *fit;
    od_rtcm_solver_t *solver;
    int source; /* 0 unset, 1 NMEA, 2 RTCM */
    char nmea[OD_MAX_SENTENCE + 1];
    size_t nmea_length;
    uint8_t rtcm[1029];
    size_t rtcm_length;
    int64_t nmea_reference, latest;
};

static int64_t od_days_from_civil(int year, unsigned month, unsigned day)
{
    int era;
    unsigned yoe, doy, doe;
    year -= month <= 2;
    era = (year >= 0 ? year : year - 399) / 400;
    yoe = (unsigned)(year - era * 400);
    doy = (153 * (month > 2 ? month - 3 : month + 9) + 2) / 5 + day - 1;
    doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return (int64_t)era * 146097 + doe - 719468;
}

static od_date_t od_civil_from_days(int64_t days)
{
    int64_t era;
    unsigned doe, yoe, doy, mp, day, month;
    int year;
    od_date_t result;
    days += 719468;
    era = (days >= 0 ? days : days - 146096) / 146097;
    doe = (unsigned)(days - era * 146097);
    yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    year = (int)yoe + (int)era * 400;
    doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    mp = (5 * doy + 2) / 153;
    day = doy - (153 * mp + 2) / 5 + 1;
    month = mp < 10 ? mp + 3 : mp - 9;
    year += month <= 2;
    result.year = year; result.month = (int)month; result.day = (int)day;
    return result;
}

static od_fit_time_t od_make_time(int64_t milliseconds)
{
    int64_t days = milliseconds / OD_DAY_MS;
    int64_t tod = milliseconds % OD_DAY_MS;
    od_date_t date;
    od_fit_time_t result = {0};
    if (tod < 0) { tod += OD_DAY_MS; --days; }
    date = od_civil_from_days(days);
    result.year = date.year; result.month = date.month; result.day = date.day;
    result.hour = (int)(tod / 3600000);
    result.minute = (int)((tod % 3600000) / 60000);
    result.second = (double)(tod % 60000) / 1000.0;
    result.unix_seconds = (double)milliseconds / 1000.0;
    result.jd_utc = result.unix_seconds / 86400.0 + 2440587.5;
    return result;
}

static int od_parse_int(const char *text, int *value)
{
    char *end;
    long parsed;
    if (!text || !*text) return 0;
    errno = 0;
    parsed = strtol(text, &end, 10);
    if (end == text || *end || errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX) return 0;
    *value = (int)parsed;
    return 1;
}

static int od_parse_double(const char *text, double *value)
{
    char *end;
    double parsed;
    if (!text || !*text) return 0;
    errno = 0;
    parsed = strtod(text, &end);
    if (end == text || *end || errno == ERANGE || !isfinite(parsed)) return 0;
    *value = parsed;
    return 1;
}

static int od_two_digits(const char *text)
{
    if (text[0] < '0' || text[0] > '9' || text[1] < '0' || text[1] > '9') return -1;
    return (text[0]-'0')*10 + text[1]-'0';
}

static int od_parse_tod(const char *text, int64_t *tod)
{
    int h, m;
    double s;
    if (strlen(text) < 6) return 0;
    h = od_two_digits(text); m = od_two_digits(text+2);
    if (h < 0 || h > 23 || m < 0 || m > 59 || !od_parse_double(text+4, &s) || s < 0 || s >= 60) return 0;
    *tod = (int64_t)h*3600000 + (int64_t)m*60000 + (int64_t)llround(s*1000.0);
    return *tod < OD_DAY_MS;
}

static int od_parse_date(const char *text, od_date_t *date)
{
    int d, m, y;
    od_date_t check;
    if (strlen(text) != 6) return 0;
    d = od_two_digits(text); m = od_two_digits(text+2); y = od_two_digits(text+4);
    if (d < 1 || d > 31 || m < 1 || m > 12 || y < 0) return 0;
    y += y >= 80 ? 1900 : 2000;
    check = od_civil_from_days(od_days_from_civil(y, (unsigned)m, (unsigned)d));
    if (check.year != y || check.month != m || check.day != d) return 0;
    *date = check;
    return 1;
}

static int od_parse_coordinate(const char *text, const char *hemisphere, int latitude, double *value)
{
    double raw, whole, minutes, result;
    char h = hemisphere[0];
    if (!h || hemisphere[1] || (latitude ? (h != 'N' && h != 'S') : (h != 'E' && h != 'W')) ||
        !od_parse_double(text, &raw) || raw < 0) return 0;
    whole = floor(raw / 100.0); minutes = raw - whole * 100.0;
    if (minutes < 0 || minutes >= 60) return 0;
    result = whole + minutes / 60.0;
    if (result > (latitude ? 90.0 : 180.0)) return 0;
    *value = (h == 'S' || h == 'W') ? -result : result;
    return 1;
}

static od_fit_vec3_t od_lla_to_ecef(double latitude, double longitude, double height)
{
    const double pi = 3.141592653589793238462643383279502884;
    const double f = 1.0 / 298.257223563, e2 = f * (2.0-f);
    double lat = latitude*pi/180.0, lon = longitude*pi/180.0;
    double sl = sin(lat), cl = cos(lat);
    double n = 6378137.0 / sqrt(1.0-e2*sl*sl);
    od_fit_vec3_t p = {(n+height)*cl*cos(lon), (n+height)*cl*sin(lon), (n*(1.0-e2)+height)*sl};
    return p;
}

static void od_push_position(od_context_t *ctx, int64_t time, od_fit_vec3_t position, od_feed_info_t *info)
{
    od_fit_observation_t observation;
    if (time < 0 || time > OD_MAX_TIME || time <= ctx->latest) {
        ++info->rejected_records; return;
    }
    observation.time_utc = od_make_time(time);
    observation.r_ecef_m = position;
    if (od_fit_context_push(ctx->fit, &observation) != OD_FIT_OK) {
        ++info->rejected_records; return;
    }
    ctx->latest = time;
    ++info->accepted_observations;
}

static int od_hex(char c)
{
    if (c >= '0' && c <= '9') return c-'0';
    if (c >= 'a' && c <= 'f') return c-'a'+10;
    if (c >= 'A' && c <= 'F') return c-'A'+10;
    return -1;
}

static void od_process_sentence(od_context_t *ctx, od_feed_info_t *info)
{
    char *fields[64], *p, *star;
    size_t n = 1, length = ctx->nmea_length, i;
    int quality;
    int64_t tod, time, base, best, distance;
    double latitude, longitude, height, geoid = 0;
    uint32_t before;
    od_date_t date;
    ctx->nmea[length] = '\0';
    while (length && isspace((unsigned char)ctx->nmea[length-1])) ctx->nmea[--length] = '\0';
    star = strchr(ctx->nmea, '*');
    if (star) {
        uint8_t checksum = 0;
        int high, low;
        if (strlen(star) != 3 || (high = od_hex(star[1])) < 0 || (low = od_hex(star[2])) < 0) goto reject;
        for (p = ctx->nmea+1; p < star; ++p) checksum ^= (uint8_t)*p;
        if (checksum != (uint8_t)((high << 4) | low)) goto reject;
        *star = '\0';
    }
    fields[0] = ctx->nmea;
    for (p = ctx->nmea; *p; ++p) if (*p == ',') {
        if (n == sizeof(fields)/sizeof(fields[0])) goto reject;
        *p = '\0'; fields[n++] = p+1;
    }
    if (strlen(fields[0]) != 6) goto reject;
    if (!strcmp(fields[0]+3, "RMC")) {
        if (n < 10 || strcmp(fields[2], "A") || !od_parse_date(fields[9], &date) ||
            !od_parse_tod(fields[1], &tod)) goto reject;
        time = od_days_from_civil(date.year, (unsigned)date.month, (unsigned)date.day)*OD_DAY_MS + tod;
        if (time < ctx->latest) goto reject;
        ctx->nmea_reference = time;
        return;
    }
    if (strcmp(fields[0]+3, "GGA")) return;
    if (n < 12 || ctx->nmea_reference < 0 || !od_parse_int(fields[6], &quality) || quality <= 0 ||
        !od_parse_tod(fields[1], &tod) || !od_parse_coordinate(fields[2], fields[3], 1, &latitude) ||
        !od_parse_coordinate(fields[4], fields[5], 0, &longitude) || !od_parse_double(fields[9], &height) ||
        strcmp(fields[10], "M") || (fields[11][0] && (!od_parse_double(fields[11], &geoid) ||
        n < 13 || strcmp(fields[12], "M")))) goto reject;
    base = (ctx->nmea_reference / OD_DAY_MS)*OD_DAY_MS + tod;
    best = base-OD_DAY_MS; distance = llabs(best-ctx->nmea_reference);
    for (i = 0; i < 2; ++i) {
        time = base + (int64_t)i*OD_DAY_MS;
        if (llabs(time-ctx->nmea_reference) < distance) {
            best = time; distance = llabs(time-ctx->nmea_reference);
        }
    }
    before = info->accepted_observations;
    od_push_position(ctx, best, od_lla_to_ecef(latitude, longitude, height+geoid), info);
    if (before != info->accepted_observations) ctx->nmea_reference = best;
    return;
reject:
    ++info->rejected_records;
}

static uint32_t od_unsigned_bits(const uint8_t *frame, size_t start, size_t length)
{
    uint32_t result = 0;
    size_t i;
    for (i = 0; i < length; ++i)
        result = (result << 1) | ((frame[(start+i)/8] >> (7-(start+i)%8)) & 1);
    return result;
}

static void od_rtcm_ephemeris_time(od_context_t *ctx, const uint8_t *frame, size_t size,
                                   unsigned type, int correct)
{
    rtcm_t *decoder = &ctx->solver->decoder;
    int week;
    double toe, toc;
    eph_t *ephemeris;
    if (type == 1019 && size*8 >= 328) {
        week = (int)od_unsigned_bits(frame,42,10) + (int)ctx->config.gps_week_rollover;
        toe = (double)od_unsigned_bits(frame,312,16)*16.0;
        toc = (double)od_unsigned_bits(frame,80,16)*16.0;
        if (!correct) decoder->time = gpst2time(week,toe);
    } else if (type == 1042 && size*8 >= 340) {
        week = (int)od_unsigned_bits(frame,42,13);
        toe = (double)od_unsigned_bits(frame,323,17)*8.0;
        toc = (double)od_unsigned_bits(frame,78,17)*8.0;
        if (!correct) decoder->time = bdt2gpst(bdt2time(week,toe));
    } else return;
    if (!correct) { ctx->solver->has_time = 1; return; }
    if (decoder->ephsat <= 0 || decoder->ephsat > MAXSAT) return;
    ephemeris = &decoder->nav.eph[decoder->ephsat-1];
    ephemeris->week = week;
    ephemeris->toe = type == 1019 ? gpst2time(week,ephemeris->toes) : bdt2gpst(bdt2time(week,ephemeris->toes));
    ephemeris->toc = type == 1019 ? gpst2time(week,toc) : bdt2gpst(bdt2time(week,toc));
    ephemeris->ttr = decoder->time;
}

static void od_decode_frame(od_context_t *ctx, const uint8_t *frame, size_t size, od_feed_info_t *info)
{
    unsigned type = ((unsigned)frame[3] << 4) | (frame[4] >> 4);
    int ephemeris_decoded = 0;
    size_t i;
    rtcm_t *decoder = &ctx->solver->decoder;
    od_rtcm_ephemeris_time(ctx,frame,size,type,0);
    if (type >= 1071 && type <= 1137 && !ctx->solver->has_time) return;
    for (i = 0; i < size; ++i) {
        int result = input_rtcm3(decoder,frame[i]);
        sol_t solution = {0};
        ssat_t satellite_status[MAXSAT] = {{0}};
        double azel[MAXOBS*2] = {0};
        char message[256] = {0};
        gtime_t utc;
        od_fit_vec3_t position;
        if (result == -1) { ++info->rejected_records; return; }
        if (result == 2) ephemeris_decoded = 1;
        if (result != 1) continue;
        if (!pntpos(decoder->obs.data,decoder->obs.n,&decoder->nav,&ctx->solver->options,
                    &solution,azel,satellite_status,message)) continue;
        if (solution.stat == SOLQ_NONE || solution.ns < 4 || !isfinite(solution.rr[0]) ||
            !isfinite(solution.rr[1]) || !isfinite(solution.rr[2])) continue;
        utc = gpst2utc(solution.time);
        position.x = solution.rr[0]; position.y = solution.rr[1]; position.z = solution.rr[2];
        od_push_position(ctx,(int64_t)utc.time*1000+(int64_t)llround(utc.sec*1000),position,info);
        ctx->solver->has_time = 1;
    }
    if (ephemeris_decoded) od_rtcm_ephemeris_time(ctx,frame,size,type,1);
}

od_config_t od_default_config(void)
{
    od_config_t config = {3,2,-1,2048};
    return config;
}

od_status_t od_create(const od_config_t *config, od_context_t **out_ctx)
{
    od_config_t cfg = config ? *config : od_default_config();
    od_context_t *ctx;
    od_fit_options_t options;
    if (!out_ctx) return OD_ERROR_ARGUMENT;
    *out_ctx = NULL;
    if (cfg.observation_capacity < 2 || cfg.observation_capacity > 65536 ||
        cfg.fit_degree < 1 || cfg.fit_degree > OD_FIT_MAX_DEGREE || cfg.fit_degree >= cfg.observation_capacity ||
        cfg.nmea_reference_utc_ms < -1 || cfg.nmea_reference_utc_ms > OD_MAX_TIME ||
        cfg.gps_week_rollover > 8192 || cfg.gps_week_rollover % 1024) return OD_ERROR_ARGUMENT;
    ctx = (od_context_t *)calloc(1,sizeof(*ctx));
    if (!ctx) return OD_ERROR_MEMORY;
    ctx->config = cfg; ctx->nmea_reference = cfg.nmea_reference_utc_ms; ctx->latest = -1;
    ctx->observations = (od_fit_observation_t *)calloc(cfg.observation_capacity,sizeof(*ctx->observations));
    options.degree = (int)cfg.fit_degree;
    if (!ctx->observations || od_fit_context_create(&ctx->fit,ctx->observations,cfg.observation_capacity,&options) != OD_FIT_OK) {
        free(ctx->observations); free(ctx); return OD_ERROR_MEMORY;
    }
    *out_ctx = ctx;
    return OD_OK;
}

od_status_t od_reset(od_context_t *ctx)
{
    if (!ctx) return OD_ERROR_ARGUMENT;
    if (ctx->solver) { free_rtcm(&ctx->solver->decoder); free(ctx->solver); ctx->solver = NULL; }
    od_fit_context_reset(ctx->fit);
    ctx->source = 0; ctx->nmea_length = ctx->rtcm_length = 0;
    ctx->nmea_reference = ctx->config.nmea_reference_utc_ms; ctx->latest = -1;
    return OD_OK;
}

void od_destroy(od_context_t *ctx)
{
    if (!ctx) return;
    od_reset(ctx);
    od_fit_context_destroy(ctx->fit);
    free(ctx->observations); free(ctx);
}

static od_status_t od_validate_feed(od_context_t *ctx, const void *data, size_t length,
                                     int source, od_feed_info_t *info)
{
    memset(info,0,sizeof(*info));
    if (!ctx || (!data && length) || length > OD_MAX_CHUNK) return OD_ERROR_ARGUMENT;
    info->buffered_observations = (uint32_t)od_fit_context_count(ctx->fit);
    if (!length) return OD_OK;
    if (ctx->source && ctx->source != source) return OD_ERROR_SOURCE;
    ctx->source = source;
    return OD_OK;
}

od_status_t od_feed_nmea(od_context_t *ctx, const char *data, size_t length, od_feed_info_t *out_info)
{
    od_feed_info_t info;
    od_status_t status = od_validate_feed(ctx,data,length,1,&info);
    size_t i;
    if (status != OD_OK || !length) { if (out_info) *out_info = info; return status; }
    for (i = 0; i < length; ++i) {
        unsigned char c = (unsigned char)data[i];
        char *star;
        if (c == '$') {
            if (ctx->nmea_length) ++info.rejected_records;
            ctx->nmea[0] = '$'; ctx->nmea[1] = '\0'; ctx->nmea_length = 1;
            continue;
        }
        if (!ctx->nmea_length) continue;
        if (c == '\r' || c == '\n') {
            od_process_sentence(ctx,&info); ctx->nmea_length = 0; continue;
        }
        if (ctx->nmea_length >= OD_MAX_SENTENCE || c < 32 || c > 126) {
            ++info.rejected_records; ctx->nmea_length = 0; continue;
        }
        ctx->nmea[ctx->nmea_length++] = (char)c; ctx->nmea[ctx->nmea_length] = '\0';
        star = strchr(ctx->nmea,'*');
        if (star && ctx->nmea_length == (size_t)(star-ctx->nmea)+3) {
            od_process_sentence(ctx,&info); ctx->nmea_length = 0;
        }
    }
    info.buffered_observations = (uint32_t)od_fit_context_count(ctx->fit);
    if (out_info) *out_info = info;
    return OD_OK;
}

static void od_remove_rtcm_prefix(od_context_t *ctx, size_t count)
{
    ctx->rtcm_length -= count;
    memmove(ctx->rtcm,ctx->rtcm+count,ctx->rtcm_length);
}

od_status_t od_feed_rtcm(od_context_t *ctx, const uint8_t *data, size_t length, od_feed_info_t *out_info)
{
    od_feed_info_t info;
    od_status_t status = od_validate_feed(ctx,data,length,2,&info);
    size_t i;
    if (status != OD_OK || !length) { if (out_info) *out_info = info; return status; }
    if (!ctx->solver) {
        od_rtcm_solver_t *s = (od_rtcm_solver_t *)calloc(1,sizeof(*s));
        if (!s || !init_rtcm(&s->decoder)) {
            free(s); if (out_info) *out_info = info; return OD_ERROR_MEMORY;
        }
        s->options = prcopt_default;
        s->options.mode = PMODE_SINGLE; s->options.navsys = SYS_GPS | SYS_CMP; s->options.nf = 1;
        s->options.elmin = 10.0*D2R; s->options.ionoopt = IONOOPT_BRDC;
        s->options.tropopt = TROPOPT_SAAS; s->options.sateph = EPHOPT_BRDC;
        ctx->solver = s;
    }
    for (i = 0; i < length; ++i) {
        ctx->rtcm[ctx->rtcm_length++] = data[i];
        for (;;) {
            size_t start = 0, payload, frame_size;
            uint32_t expected;
            while (start < ctx->rtcm_length && ctx->rtcm[start] != 0xD3) ++start;
            od_remove_rtcm_prefix(ctx,start);
            if (ctx->rtcm_length < 3) break;
            if (ctx->rtcm[1] & 0xFC) {
                ++info.rejected_records; od_remove_rtcm_prefix(ctx,1); continue;
            }
            payload = ((ctx->rtcm[1] & 3U) << 8) | ctx->rtcm[2];
            frame_size = payload+6;
            if (ctx->rtcm_length < frame_size) break;
            expected = ((uint32_t)ctx->rtcm[frame_size-3] << 16) |
                       ((uint32_t)ctx->rtcm[frame_size-2] << 8) | ctx->rtcm[frame_size-1];
            if (rtk_crc24q(ctx->rtcm,(int)frame_size-3) != expected) {
                ++info.rejected_records; od_remove_rtcm_prefix(ctx,1); continue;
            }
            if (payload < 20) ++info.rejected_records;
            else od_decode_frame(ctx,ctx->rtcm,frame_size,&info);
            od_remove_rtcm_prefix(ctx,frame_size);
        }
    }
    info.buffered_observations = (uint32_t)od_fit_context_count(ctx->fit);
    if (out_info) *out_info = info;
    return OD_OK;
}

od_status_t od_get_orbit(od_context_t *ctx, od_j2000_state_t *out_state)
{
    od_fit_time_t time;
    od_fit_state_t state;
    od_j2000_state_t orbit;
    int i;
    if (!ctx || !out_state) return OD_ERROR_ARGUMENT;
    if (od_fit_context_count(ctx->fit) < ctx->config.fit_degree+1) return OD_NOT_READY;
    time = od_make_time(ctx->latest);
    if (od_fit_context_query_state(ctx->fit,&time,&state) != OD_FIT_OK) return OD_ERROR_FIT;
    orbit.timestamp_utc_ms = ctx->latest;
    orbit.position_m[0] = state.r_j2000_m.x; orbit.position_m[1] = state.r_j2000_m.y; orbit.position_m[2] = state.r_j2000_m.z;
    orbit.velocity_mps[0] = state.v_j2000_mps.x; orbit.velocity_mps[1] = state.v_j2000_mps.y; orbit.velocity_mps[2] = state.v_j2000_mps.z;
    for (i = 0; i < 3; ++i)
        if (!isfinite(orbit.position_m[i]) || !isfinite(orbit.velocity_mps[i])) return OD_ERROR_FIT;
    *out_state = orbit;
    return OD_OK;
}

const char *od_status_string(od_status_t status)
{
    switch (status) {
    case OD_OK: return "ok";
    case OD_NOT_READY: return "not enough observations";
    case OD_ERROR_ARGUMENT: return "invalid argument or configuration";
    case OD_ERROR_MEMORY: return "allocation failed";
    case OD_ERROR_SOURCE: return "reset before changing input source";
    case OD_ERROR_FIT: return "orbit fit failed";
    case OD_ERROR_INTERNAL: return "internal error; reset before retrying";
    default: return "unknown status";
    }
}
