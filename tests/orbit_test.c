/* Pure C tests; include the single implementation for private numeric checks. */
#include "../src/orbit_determination.c"
#include <inttypes.h>

#define CHECK(x) do { if (!(x)) { fprintf(stderr,"line %d: %s\n",__LINE__,#x); exit(1); } } while (0)

static od_context_t *test_create(uint32_t capacity, int64_t reference)
{
    od_config_t config = od_default_config();
    od_context_t *ctx = NULL;
    config.observation_capacity = capacity;
    config.fit_degree = 2; /* Existing streaming regressions use quadratic fits. */
    config.nmea_reference_utc_ms = reference;
    CHECK(od_create(&config,&ctx) == OD_OK);
    return ctx;
}

static void test_checked(const char *body, char *out, size_t capacity)
{
    uint8_t checksum = 0;
    const char *p;
    int n;
    for (p=body; *p; ++p) checksum ^= (uint8_t)*p;
    n=snprintf(out,capacity,"$%s*%02X\r\n",body,(unsigned)checksum);
    CHECK(n > 0 && (size_t)n < capacity);
}

static void test_rmc(char *out)
{
    test_checked("GNRMC,040403.00,A,1724.2594222,N,15304.2955483,W,0.0,0.0,230626,,,A",out,256);
}

static void test_gga_time(const char *time, int point, char *out)
{
    char body[256];
    snprintf(body,sizeof(body),"GNGGA,%s,1724.2594222,N,%.7f,W,1,14,0.69,271441.347,M,0.989,M,,",
             time,15304.2955483+point*0.02);
    test_checked(body,out,256);
}

static void test_gga(int point, char *out)
{
    char time[32];
    snprintf(time,sizeof(time),"0404%02d.00",3+point);
    test_gga_time(time,point,out);
}

static od_feed_info_t test_send(od_context_t *ctx, const char *data, size_t chunk)
{
    od_feed_info_t total={0};
    size_t length=strlen(data), offset;
    if (!chunk) chunk=length;
    for (offset=0; offset<length; offset+=chunk) {
        od_feed_info_t info;
        size_t count=length-offset < chunk ? length-offset : chunk;
        CHECK(od_feed_nmea(ctx,data+offset,count,&info) == OD_OK);
        total.accepted_observations+=info.accepted_observations;
        total.rejected_records+=info.rejected_records;
        total.buffered_observations=info.buffered_observations;
    }
    return total;
}

static void test_same(const od_j2000_state_t *a, const od_j2000_state_t *b)
{
    int i;
    CHECK(a->timestamp_utc_ms==b->timestamp_utc_ms);
    for (i=0;i<3;++i) {
        CHECK(fabs(a->position_m[i]-b->position_m[i])<1e-6);
        CHECK(fabs(a->velocity_mps[i]-b->velocity_mps[i])<1e-6);
    }
}

static void test_default_config(void)
{
    od_config_t cfg = od_default_config();
    od_context_t *implicit = NULL, *explicit = NULL, *tail = NULL;
    od_j2000_state_t a, b;
    char data[256];
    int i;
    CHECK(cfg.observation_capacity == 30 && cfg.fit_degree == 3);
    CHECK(od_create(NULL, &implicit) == OD_OK);
    CHECK(od_create(&cfg, &explicit) == OD_OK);
    CHECK(od_create(&cfg, &tail) == OD_OK);
    test_rmc(data);
    test_send(implicit, data, 1); test_send(explicit, data, 0); test_send(tail, data, 0);
    for (i = 0; i < 40; ++i) {
        od_feed_info_t info;
        test_gga(i, data);
        info = test_send(implicit, data, 1);
        CHECK(info.accepted_observations == 1 && info.rejected_records == 0);
        CHECK(info.buffered_observations == (uint32_t)(i < 30 ? i + 1 : 30));
        test_send(explicit, data, 0);
        CHECK(od_get_orbit(implicit, &a) == (i < 3 ? OD_NOT_READY : OD_OK));
        CHECK(od_get_orbit(explicit, &b) == (i < 3 ? OD_NOT_READY : OD_OK));
        if (i >= 3) test_same(&a, &b);
        if (i >= 10) test_send(tail, data, 0);
    }
    CHECK(od_get_orbit(tail, &b) == OD_OK); test_same(&a, &b);
    CHECK(od_reset(implicit) == OD_OK);
    test_rmc(data); test_send(implicit, data, 0);
    for (i = 0; i < 4; ++i) {
        test_gga(i, data); test_send(implicit, data, 0);
        CHECK(od_get_orbit(implicit, &a) == (i < 3 ? OD_NOT_READY : OD_OK));
    }
    od_destroy(implicit); od_destroy(explicit); od_destroy(tail);
}

static void test_api(void)
{
    od_context_t *raw=NULL, *ctx, *midnight;
    od_config_t cfg=od_default_config();
    od_j2000_state_t orbit, unchanged;
    char data[256], rmc[256], long_data[3400];
    int i, c;
    const uint32_t capacities[]={3,5,10};
    const uint8_t frame[]={0xD3,0,4,0x43,0x50,0,0,0x44,0xFE,0x2E};
    uint32_t rejected=0;
    test_default_config();
    cfg.observation_capacity=2;
    CHECK(od_create(&cfg,&raw)==OD_ERROR_ARGUMENT && !raw);
    cfg.fit_degree=1;
    CHECK(od_create(&cfg,&raw)==OD_OK);
    od_destroy(raw);
    cfg.gps_week_rollover=1;
    CHECK(od_create(&cfg,&raw)==OD_ERROR_ARGUMENT);
    CHECK(od_get_orbit(NULL,NULL)==OD_ERROR_ARGUMENT);
    CHECK(od_reset(NULL)==OD_ERROR_ARGUMENT);
    od_destroy(NULL);
    ctx=test_create(3,-1);
    memset(&orbit,0x5A,sizeof(orbit)); memcpy(&unchanged,&orbit,sizeof(orbit));
    CHECK(od_get_orbit(ctx,&orbit)==OD_NOT_READY);
    CHECK(!memcmp(&orbit,&unchanged,sizeof(orbit)));
    test_gga(0,data); CHECK(test_send(ctx,data,0).rejected_records==1);
    test_rmc(rmc); CHECK(test_send(ctx,rmc,1).accepted_observations==0);
    for (i=0;i<3;++i) {
        od_feed_info_t info;
        test_gga(i,data); info=test_send(ctx,data,1);
        CHECK(info.accepted_observations==1 && !info.rejected_records);
        CHECK(info.buffered_observations==(uint32_t)(i+1));
        CHECK(od_get_orbit(ctx,&orbit)==(i<2 ? OD_NOT_READY : OD_OK));
    }
    CHECK(orbit.timestamp_utc_ms==INT64_C(1782187445000));
    unchanged=orbit;
    test_gga(2,data); CHECK(test_send(ctx,data,0).rejected_records==1);
    test_gga(1,data); CHECK(test_send(ctx,data,0).rejected_records==1);
    CHECK(od_get_orbit(ctx,&orbit)==OD_OK); test_same(&orbit,&unchanged);
    test_gga(3,data); data[10]=data[10]=='0' ? '1' : '0';
    CHECK(test_send(ctx,data,0).rejected_records==1);
    memset(long_data,'X',3001); long_data[0]='$'; long_data[3001]='\n'; long_data[3002]=0;
    test_gga(3,data); strcat(long_data,data);
    CHECK(test_send(ctx,long_data,7).accepted_observations==1);
    CHECK(od_feed_rtcm(ctx,frame,1,NULL)==OD_ERROR_SOURCE);
    CHECK(od_feed_nmea(ctx,NULL,1,NULL)==OD_ERROR_ARGUMENT);
    CHECK(od_feed_nmea(ctx,NULL,0,NULL)==OD_OK);
    for (c=0;c<3;++c) {
        uint32_t capacity=capacities[c];
        od_context_t *full=test_create(capacity,-1), *tail=test_create(capacity,-1);
        test_send(full,rmc,0); test_send(tail,rmc,0);
        for (i=0;i<40;++i) {
            od_feed_info_t info;
            test_gga(i,data); info=test_send(full,data,0);
            CHECK(info.buffered_observations==((uint32_t)(i+1)<capacity ? (uint32_t)(i+1) : capacity));
            CHECK(od_get_orbit(full,&orbit)==(i<2 ? OD_NOT_READY : OD_OK));
            if (i>=40-(int)capacity) test_send(tail,data,11);
        }
        CHECK(od_get_orbit(tail,&unchanged)==OD_OK); test_same(&orbit,&unchanged);
        od_destroy(full); od_destroy(tail);
    }
    midnight=test_create(3,INT64_C(1782259199000));
    test_gga_time("235959.00",0,data); CHECK(test_send(midnight,data,0).accepted_observations==1);
    test_gga_time("000000.00",1,data); CHECK(test_send(midnight,data,0).accepted_observations==1);
    test_gga_time("000001.00",2,data); CHECK(test_send(midnight,data,0).accepted_observations==1);
    CHECK(od_get_orbit(midnight,&orbit)==OD_OK && orbit.timestamp_utc_ms==INT64_C(1782259201000));
    od_destroy(midnight);
    CHECK(od_reset(ctx)==OD_OK);
    for (i=0;i<(int)sizeof(frame);++i) {
        od_feed_info_t info;
        CHECK(od_feed_rtcm(ctx,frame+i,1,&info)==OD_OK); rejected+=info.rejected_records;
    }
    CHECK(rejected==1 && od_get_orbit(ctx,&orbit)==OD_NOT_READY);
    CHECK(od_reset(ctx)==OD_OK);
    test_send(ctx,rmc,0);
    for(i=0;i<3;++i) { test_gga(i,data); CHECK(test_send(ctx,data,0).accepted_observations==1); }
    CHECK(od_get_orbit(ctx,&orbit)==OD_OK);
    od_destroy(ctx);
    puts("pure C API, streaming, time and ring checks passed");
}

static od_j2000_state_t test_replay(od_context_t *ctx, const uint8_t *data, size_t length, size_t chunk, uint32_t *total)
{
    size_t i;
    od_j2000_state_t result;
    *total=0;
    for(i=0;i<length;i+=chunk) {
        od_feed_info_t info;
        size_t count=length-i<chunk ? length-i : chunk;
        CHECK(od_feed_rtcm(ctx,data+i,count,&info)==OD_OK);
        *total+=info.accepted_observations; CHECK(info.buffered_observations<=3);
    }
    CHECK(*total>=3 && od_get_orbit(ctx,&result)==OD_OK);
    return result;
}

static void test_capture(const char *path)
{
    FILE *file=fopen(path,"rb");
    long length;
    uint8_t *data;
    od_context_t *ctx=test_create(3,-1);
    od_j2000_state_t reference;
    uint32_t total;
    size_t chunks[]={1,1024,1048576}, i;
    CHECK(file && fseek(file,0,SEEK_END)==0);
    length=ftell(file); CHECK(length>0 && fseek(file,0,SEEK_SET)==0);
    data=(uint8_t *)malloc((size_t)length); CHECK(data);
    CHECK(fread(data,1,(size_t)length,file)==(size_t)length); fclose(file);
    reference=test_replay(ctx,data,(size_t)length,257,&total);
    for(i=0;i<3;++i) {
        uint32_t count;
        od_j2000_state_t result;
        CHECK(od_reset(ctx)==OD_OK);
        result=test_replay(ctx,data,(size_t)length,chunks[i],&count);
        CHECK(count==total); test_same(&result,&reference);
    }
    free(data); od_destroy(ctx);
    printf("RTCM capture: %u positions; fragmentation/reset invariant\n",total);
}

/* Numeric core regression follows. */
static int test_math(void)
{
    od_fit_observation_t buffer[3];
    od_fit_context_t *ctx = NULL;
    od_fit_options_t options = {2};
    const double epoch = 1782187443.0;
    int i, j;
    CHECK(od_fit_context_create(&ctx, buffer, 3, &options) == OD_FIT_OK);
    for (i = 0; i < 20; ++i) {
        const double t = i * 0.7; /* Supports non-integer-second epochs. */
        const double r[3] = {7000000 + 20*t + 0.5*t*t, 7000*t - 2*t*t, 1000 - 3*t + 4*t*t};
        const double v[3] = {20+t, 7000-4*t, -3+8*t};
        double matrix[3][3];
        od_fit_observation_t obs;
        od_fit_state_t state;
        od_fit_time_t future;
        memset(&obs, 0, sizeof(obs));
        obs.time_utc.unix_seconds = epoch + t;
        obs.time_utc.jd_utc = (epoch + t) / 86400 + 2440587.5;
        od_fit_itrs_to_j2000_matrix(&obs.time_utc, matrix);
        obs.r_ecef_m.x = matrix[0][0]*r[0] + matrix[1][0]*r[1] + matrix[2][0]*r[2];
        obs.r_ecef_m.y = matrix[0][1]*r[0] + matrix[1][1]*r[1] + matrix[2][1]*r[2];
        obs.r_ecef_m.z = matrix[0][2]*r[0] + matrix[1][2]*r[1] + matrix[2][2]*r[2];
        CHECK(od_fit_context_push(ctx, &obs) == OD_FIT_OK);
        CHECK(od_fit_context_count(ctx) == (size_t)(i < 3 ? i+1 : 3));
        if (i < 2) continue;
        CHECK(od_fit_context_query_state(ctx, &obs.time_utc, &state) == OD_FIT_OK);
        {
            double actual_r[3] = {state.r_j2000_m.x,state.r_j2000_m.y,state.r_j2000_m.z};
            double actual_v[3] = {state.v_j2000_mps.x,state.v_j2000_mps.y,state.v_j2000_mps.z};
            for (j = 0; j < 3; ++j) {
                CHECK(fabs(actual_r[j]-r[j]) < 1e-5);
                CHECK(fabs(actual_v[j]-v[j]) < 0.01);
            }
        }
        future = obs.time_utc;
        future.unix_seconds += 0.001;
        CHECK(od_fit_context_query_state(ctx, &future, &state) == OD_FIT_ERR_RANGE);
        future.unix_seconds = epoch - 1;
        CHECK(od_fit_context_query_state(ctx, &future, &state) == OD_FIT_ERR_RANGE);
        CHECK(od_fit_context_push(ctx, &obs) == OD_FIT_ERR_RANGE);
    }
    od_fit_context_destroy(ctx);
    puts("numeric fitting, ring wrap and no-extrapolation checks passed");
    return 0;
}

int main(int argc, char **argv)
{
    if (argc==2 && !strcmp(argv[1],"api")) { test_api(); return 0; }
    if (argc==2 && !strcmp(argv[1],"core")) return test_math();
    if (argc==3 && !strcmp(argv[1],"rtcm")) { test_capture(argv[2]); return 0; }
    fprintf(stderr,"usage: orbit_test api|core|rtcm FILE\n");
    return 2;
}
