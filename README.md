# calgnss_orbit_mini

`calgnss_orbit_mini` 是一个面向嵌入式/跨平台集成的 C 语言小型轨道拟合库。库接收按时间递增的 ECEF/ITRS 位置观测数据，在内部维护一段可配置长度的环形缓存，并根据用户查询时间输出 J2000 惯性坐标系下的位置和速度状态量。查询时间位于观测范围内时执行插值拟合，位于最新观测之后时支持短期未来轨道外推。

## 设计目标

- 开发语言：C99。
- 依赖范围：仅使用 C 标准库头文件和 `math.h`。
- 不依赖 Windows API、POSIX API、文件系统、线程、动态加载或系统时间函数。
- 观测缓存由调用方提供，缓存长度 `N` 可配置。
- 面向 Windows、Linux、RTOS 等平台移植。
- 目标资源约束：RAM 尽量控制在 128 KB 以内，CPU 消耗尽量低。

说明：当前实现的 `cg_context_create()` 会对上下文结构体执行一次 `malloc()`，`cg_context_destroy()` 会执行一次 `free()`；观测环形缓存由调用方静态或外部管理。如果目标 RTOS 严格禁止堆内存，可将上下文结构体改为由调用方静态分配，算法主体不依赖文件、线程或系统时间。

## 输入

每组观测数据使用 `cg_observation_t` 表示：

```c
typedef struct cg_observation_t {
    cg_time_t time_utc;
    cg_vec3_t r_ecef_m;
} cg_observation_t;
```

字段含义：

| 字段 | 单位 | 说明 |
| --- | --- | --- |
| `time_utc` | UTC | 观测时间。需要同时填写日历时间、`jd_utc` 和 `unix_seconds`。 |
| `r_ecef_m.x/y/z` | m | 地固系 ECEF/ITRS 位置。 |

输入要求：

- 数据必须按时间严格递增输入，推荐时间间隔为 1 s。
- 允许存在不定时长的数据缺失，但查询精度会随缺口长度、观测质量和缓存覆盖范围下降。
- 坐标基准必须为 ECEF/ITRS，单位为 m。
- 调用方负责把外部观测数据解析为 `cg_observation_t`，本库不读文件、不访问串口、不依赖系统时间。

时间结构：

```c
typedef struct cg_time_t {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
    double jd_utc;
    double unix_seconds;
} cg_time_t;
```

库内部主要使用：

- `jd_utc`：用于 ECEF 到 J2000 的时间相关坐标转换。
- `unix_seconds`：用于环形缓存排序、插值拟合和时间范围检查。

因此，调用方必须保证日历时间、`jd_utc` 和 `unix_seconds` 一致。

## 输出

查询输出使用 `cg_state_t` 表示：

```c
typedef struct cg_state_t {
    cg_time_t time_utc;
    cg_vec3_t r_j2000_m;
    cg_vec3_t v_j2000_mps;
} cg_state_t;
```

字段含义：

| 字段 | 单位 | 说明 |
| --- | --- | --- |
| `time_utc` | UTC | 查询时间。 |
| `r_j2000_m.x/y/z` | m | J2000 惯性坐标系下的位置。 |
| `v_j2000_mps.x/y/z` | m/s | J2000 惯性坐标系下的速度。 |

输出能力：

- 查询时间位于已有观测时间范围内时，输出该时刻的 J2000 状态量。查询时间可以不是原始输入数据中的时间点，支持毫秒级时间输入和插值拟合。
- 查询时间晚于最新观测时间时，支持最长 3600 s 的短期未来轨道外推；超过最大外推时长时返回 `CG_ERR_RANGE`。

目标精度：

- 观测时间范围内插值：定位精度目标 `< 10 m`，测速精度目标 `< 0.2 m/s`。
- 最新观测后的短期外推：精度取决于观测弧长、观测质量、缺测情况、轨道环境和外推时长，建议在目标场景中用实测数据验证。

精度前提：

- 输入 ECEF 观测本身的时间和坐标误差满足目标精度要求。
- 缓存中有足够连续、覆盖合理的观测数据。
- 查询时间必须落在缓存观测时间范围内，或不超过最新观测后 3600 s。

## 处理流程

1. 调用方创建长度为 `N` 的 `cg_observation_t` 缓冲区。
2. `cg_context_create()` 绑定该缓冲区并创建上下文。
3. 按时间递增调用 `cg_context_push()` 写入观测数据。
4. 库内部以环形缓存保存最近 `N` 组观测。
5. 调用 `cg_context_query_state()` 查询任意支持时间点的状态量。
6. 查询时间在观测范围内时，使用局部 Chebyshev 最小二乘拟合计算位置和速度。
7. 查询时间晚于最新观测且不超过 3600 s 时，基于末端观测拟合最新 J2000 状态，并使用简化动力学模型进行数值积分外推。
8. 查询时间早于缓存最早观测，或超过最大未来外推时长时返回 `CG_ERR_RANGE`。

## 使用条件

### 平台条件

- 编译器需支持 C99。
- 需要支持 `double` 浮点运算。
- 需要提供 `math.h` 中的基础数学函数，例如 `sin()`、`cos()`、`floor()`、`fmod()`、`isfinite()`。
- 目标平台无需提供文件系统、线程、动态库加载、系统时间或 OS API。

### 数据条件

- 至少输入 2 组观测后才能查询。
- 观测时间必须严格递增；重复时间或倒序输入会返回 `CG_ERR_RANGE`。
- 查询时间早于缓存最早观测时间，或晚于最新观测后 3600 s 会返回 `CG_ERR_RANGE`。
- 缺测时间越长，拟合窗口内可用数据越少，插值精度越容易下降。
- 外推会优先使用最新观测前约 600 s 的历史数据估计末端状态和经验力模型；有效历史弧段越短、缺测越多，外推结果越需要谨慎使用。
- 针对 270 km 晨昏轨道，建议缓存长度 `N = 300 ~ 600`，即 1 Hz 输入下保留约 5 到 10 分钟末端观测弧，用于校准短期外推力模型。

### 缓存长度 `N`

缓存长度由调用方传入：

```c
cg_observation_t obs_buffer[N];
cg_context_create(&ctx, obs_buffer, N, &opt);
```

选择建议：

- `N >= 2`：接口最低要求，仅适合功能连通性测试。
- `N = 10`：库默认值，适合极低内存 demo，不建议作为最终精度配置。
- `N = 300 ~ 1200`：适合 1 s 输入下保存 5 到 20 分钟观测，通常是嵌入式场景更现实的起点。

典型 64 位平台上，`sizeof(cg_observation_t)` 约为 72 字节，因此观测缓存内存约为：

```text
buffer_bytes ~= N * sizeof(cg_observation_t)
```

例如 `N = 1024` 时观测缓存约 72 KB，再加上上下文、拟合缓存和函数栈临时数组，通常仍可控制在 128 KB 目标附近。不同编译器和 ABI 的结构体对齐可能不同，最终应以目标平台 `sizeof()` 结果为准。

### 选项配置

```c
typedef struct cg_options_t {
    int degree;
} cg_options_t;
```

| 选项 | 说明 |
| --- | --- |
| `degree` | Chebyshev 拟合阶数，内部上限为 `CG_MAX_DEGREE`。默认 10。 |

## 返回码

| 返回码 | 含义 |
| --- | --- |
| `CG_OK` | 成功。 |
| `CG_ERR_INVALID_ARGUMENT` | 输入指针为空、容量不足、数据量不足等非法参数。 |
| `CG_ERR_RANGE` | 时间范围不支持，例如倒序输入、查询早于缓存或查询超过最大未来外推时长。 |
| `CG_ERR_FIT` | 拟合失败，常见原因是数据量不足或几何条件差。 |
| `CG_ERR_NO_MEMORY` | 创建上下文时内存分配失败。 |
| `CG_ERR_PARSE` / `CG_ERR_IO` | 预留错误码，核心库当前不执行解析和 I/O。 |

可使用 `cg_status_string(status)` 获取简短错误描述。

## 使用 Demo

下面 demo 展示纯内存调用流程：构造时间、创建环形缓存、写入观测、查询毫秒级插值点和短期未来外推点。示例中的 ECEF 数据只用于演示 API，不代表真实轨道数据。

```c
#include "calgnss_orbit_mini.h"

#include <math.h>
#include <stdio.h>
#include <stdint.h>

#define OBS_CAPACITY 1024

static int64_t demo_days_from_civil(int y, unsigned m, unsigned d)
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

static void demo_civil_from_days(int64_t z, int *y, unsigned *m, unsigned *d)
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

static double demo_jd_from_calendar(
    int y, int m, int d, int hour, int minute, double second)
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

static cg_time_t demo_make_time(
    int y, int m, int d, int hour, int minute, double second)
{
    cg_time_t t;
    int whole_second = (int)floor(second);
    double frac = second - (double)whole_second;
    int64_t days = demo_days_from_civil(y, (unsigned)m, (unsigned)d);

    t.year = y;
    t.month = m;
    t.day = d;
    t.hour = hour;
    t.minute = minute;
    t.second = second;
    t.jd_utc = demo_jd_from_calendar(y, m, d, hour, minute, second);
    t.unix_seconds = (double)days * 86400.0 +
                     (double)hour * 3600.0 +
                     (double)minute * 60.0 +
                     (double)whole_second + frac;

    return t;
}

static cg_time_t demo_add_seconds(const cg_time_t *base, double seconds)
{
    double unix_seconds = base->unix_seconds + seconds;
    double whole_d = floor(unix_seconds);
    int64_t whole = (int64_t)whole_d;
    double frac = unix_seconds - whole_d;
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

    demo_civil_from_days(days, &y, &m, &d);
    hour = (int)(sod / 3600);
    minute = (int)((sod % 3600) / 60);

    return demo_make_time(
        y, (int)m, (int)d, hour, minute, (double)(sod % 60) + frac);
}

int main(void)
{
    static cg_observation_t obs_buffer[OBS_CAPACITY];
    cg_context_t *ctx = NULL;
    cg_options_t opt;
    cg_time_t t0;
    cg_time_t query_interp;
    cg_time_t query_future;
    cg_state_t state;
    int status;
    int i;

    opt = cg_default_options();

    status = cg_context_create(&ctx, obs_buffer, OBS_CAPACITY, &opt);
    if (status != CG_OK) {
        printf("create failed: %s\n", cg_status_string(status));
        return 1;
    }

    t0 = demo_make_time(2026, 5, 18, 7, 0, 0.0);

    for (i = 0; i < 600; ++i) {
        cg_observation_t obs;

        obs.time_utc = demo_add_seconds(&t0, (double)i);
        obs.r_ecef_m.x = 6800000.0 * cos((double)i * 0.001);
        obs.r_ecef_m.y = 6800000.0 * sin((double)i * 0.001);
        obs.r_ecef_m.z = 500000.0 * sin((double)i * 0.002);

        status = cg_context_push(ctx, &obs);
        if (status != CG_OK) {
            printf("push failed at %d: %s\n", i, cg_status_string(status));
            cg_context_destroy(ctx);
            return 1;
        }
    }

    query_interp = demo_add_seconds(&t0, 123.456);
    status = cg_context_query_state(ctx, &query_interp, &state);
    if (status == CG_OK) {
        printf("interpolation J2000 r = %.3f %.3f %.3f m\n",
               state.r_j2000_m.x, state.r_j2000_m.y, state.r_j2000_m.z);
        printf("interpolation J2000 v = %.6f %.6f %.6f m/s\n",
               state.v_j2000_mps.x, state.v_j2000_mps.y, state.v_j2000_mps.z);
    } else {
        printf("interpolation query failed: %s\n", cg_status_string(status));
    }

    query_future = demo_add_seconds(&t0, 599.0 + 120.0);
    status = cg_context_query_state(ctx, &query_future, &state);
    if (status == CG_OK) {
        printf("extrapolation J2000 r = %.3f %.3f %.3f m\n",
               state.r_j2000_m.x, state.r_j2000_m.y, state.r_j2000_m.z);
        printf("extrapolation J2000 v = %.6f %.6f %.6f m/s\n",
               state.v_j2000_mps.x, state.v_j2000_mps.y, state.v_j2000_mps.z);
    } else {
        printf("extrapolation query failed: %s\n", cg_status_string(status));
    }

    cg_context_destroy(ctx);
    return 0;
}
```

编译示例：

```sh
cmake -S . -B build
cmake --build build --config Release
```

PC 侧测试可链接动态库或静态库；嵌入式/RTOS 集成建议直接链接静态库或源码文件。

## API 快速索引

```c
cg_options_t cg_default_options(void);

int cg_context_create(
    cg_context_t **out_context,
    cg_observation_t *observation_buffer,
    size_t capacity,
    const cg_options_t *options);

void cg_context_reset(cg_context_t *context);

void cg_context_destroy(cg_context_t *context);

size_t cg_context_count(const cg_context_t *context);

size_t cg_context_capacity(const cg_context_t *context);

int cg_context_push(
    cg_context_t *context,
    const cg_observation_t *observation);

int cg_context_query_state(
    cg_context_t *context,
    const cg_time_t *query_time_utc,
    cg_state_t *out_state);

const char *cg_status_string(int status);
```

## 算法概要

- ECEF/ITRS 转换到 J2000/GCRS，内部使用紧凑的 IAU 2000B 章动/岁差和格林尼治视恒星时模型。默认 `CG_USE_EOP_TABLE=0`，按 UTC≈UT1、无极移处理，以匹配常见 Fixed/J2000 测试数据；如需启用库内置 2026 年 EOP 日表，可在编译时定义 `CG_USE_EOP_TABLE=1`。
- 观测范围内查询使用局部 Chebyshev 最小二乘位置拟合，速度由拟合多项式解析求导得到。
- 未来查询在最新观测后 3600 s 内使用末端状态估计、力模型校准和 RK4 数值积分进行短期外推；超过该范围返回 `CG_ERR_RANGE`。
- 外推动力学包含中心引力、可自适应缩放的非中心地球重力、日月三体摄动、太阳光压、大气阻力和 RTN 经验加速度。
- 大气密度由高度表、太阳活动基准、地磁活动基准和局地太阳时因子共同调制；阻力强度由末端历史弧约束，避免在无明显衰减的数据上引入虚假拖曳。
