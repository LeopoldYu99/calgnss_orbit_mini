# 定轨算法：纯 C 单文件版

对外名称统一为 `orbit_determination`，实现文件为 `src/orbit_determination.c`，头文件为 `include/orbit_determination.h`。原先的旧数学模块名称已移除，内部数学函数统一使用 `od_fit_` 前缀。

**整个库只有一个 C99 实现文件。** 数学拟合、坐标转换、EOP 数据表、NMEA 解析和所需 RTKLIB 实现均已内嵌，不通过 `#include` 拼接其他私有源码或数据文件。测试、示例和构建也只使用 C 编译器，没有 C++ 标准库依赖。`third_party/rtklib` 仅保留许可和来源说明。

功能：RTCM3 或 NMEA 输入，输出最新有效观测历元的 J2000 位置、速度及 UTC 毫秒时间戳。保留可配置的环形位置观测队列，不包含卫星自身的未来轨道外推接口、外推积分器或外推专用缓存。内嵌 RTKLIB 中的 GNSS 导航星历计算是 RTCM 定位所需功能。

## 最简集成

只需复制 `orbit_determination.c` 和 `orbit_determination.h`，使用 C 编译器即可：

```sh
gcc -std=c99 -O2 -Iinclude src/orbit_determination.c examples/od_replay.c \
    -lm -o od_replay
./od_replay rtcm /path/to/capture.rtcm3 30
./od_replay nmea /path/to/nmea.log 30
```

回放示例从文件读数据，输出最后一个状态的 CSV。省略点数参数时使用默认 30 点、三次拟合；指定点数时至少为 4，拟合阶数仍为 3。文件和通信收发由外部调用程序处理，库的输入接口接收字节。

也可构建静态库：

```sh
cmake -S . -B build-c -DCMAKE_BUILD_TYPE=Release
cmake --build build-c -j
ctest --test-dir build-c --output-on-failure
cmake --install build-c --prefix "$PWD/dist/linux-x86_64"
```

不需要联网下载依赖。CMake 项目只声明 `LANGUAGES C`。

### 不支持 `dirent.h` 的平台

源码默认 `OD_ENABLE_DIRECTORY_SCAN=0`，不包含 `dirent.h`，也不调用 `opendir/readdir/closedir`。NMEA 和 RTCM 字节流定轨不使用目录扫描，无需在硬件工程里添加此头文件或替代实现。

内嵌 RTKLIB 的旧文件路径辅助函数在此模式下只接受明确文件名，不展开 `*`/`?` 通配符。桌面平台确需使用这些辅助函数时，可在编译整个实现文件时定义 `OD_ENABLE_DIRECTORY_SCAN=1`。

### FreeRTOS / 裸机编译

源码默认采用 FreeRTOS / 裸机行为，无需定义平台宏。默认不依赖 pthread、Windows API、目录扫描、`sys/time.h`、`sys/stat.h`，不调用系统时钟、休眠、shell 或目录创建接口，也不开启 RTKLIB 文件跟踪日志。桌面编译使用相同行为；输入时间仍由 NMEA/RTCM 数据提供。

GCC 类工具链的源码编译示意（将 `gcc` 换成板卡交叉编译器，补充芯片参数）：

```sh
gcc -std=c99 -O2 -ffunction-sections -fdata-sections \
    -Iinclude -c src/orbit_determination.c -o orbit_determination.o
```

固件链接时加入 `orbit_determination.o`、数学库 `-lm` 和 `-Wl,--gc-sections`，清除不使用的内嵌 RTKLIB 文件辅助函数；无需 `-pthread`。标准 C 运行库仍需提供内存分配、字符串/数值转换、格式化和数学函数，浮点运算使用 `double`。源码不自动将 `malloc/free` 映射到 FreeRTOS 堆接口。

使用 CMake 时，在配置工具链的同时添加 `-DOD_BUILD_EXAMPLES=OFF -DOD_BUILD_TESTS=OFF`；嵌入式库不查找或链接 Threads。仓库中的 `dist/linux-x86_64/lib/libod.a` 仍是 Linux 库，硬件需用自己的工具链重新编译源码。

NMEA 时间由 RMC/GGA 或显式参考时间提供；GPS/北斗 RTCM 时间继续由输入星历与观测提供，不需要 FreeRTOS 提供当前日期。一个实例应由单个任务调用；多个 RTCM 实例也应串行解码。库不增加内部互斥锁，不改变默认 30 点、三次拟合。

## 接口

| 接口 | 作用 |
|---|---|
| `od_default_config()` | 获取完整默认配置 |
| `od_create(config, &ctx)` | 创建实例；config 为 NULL 时使用默认配置 |
| `od_feed_rtcm(ctx, data, length, &info)` | 输入原始 RTCM3 字节流 |
| `od_feed_nmea(ctx, data, length, &info)` | 输入 NMEA 文本字节流 |
| `od_get_orbit(ctx, &orbit)` | 输出最新观测历元的 J2000 位置和速度 |
| `od_reset(ctx)` | 清空输入、定位和拟合状态，保留配置 |
| `od_destroy(ctx)` | 释放实例 |
| `od_status_string(status)` | 获取静态状态说明字符串 |

初始化和输出示意：

```c
od_config_t config = od_default_config();
config.observation_capacity = 30; /* 环形队列最多保留 30 个位置观测点 */
config.fit_degree = 3;           /* 三次拟合 */

od_context_t *ctx = NULL;
od_status_t status = od_create(&config, &ctx);
/* 成功后调用 od_feed_rtcm 或 od_feed_nmea，不在同一实例混合来源。 */

od_j2000_state_t orbit;
status = od_get_orbit(ctx, &orbit);
/* OD_OK：orbit.position_m[3] 和 orbit.velocity_mps[3] 有效。
   OD_NOT_READY：有效位置观测不足。 */

od_destroy(ctx);
```

完整错误处理、数据输入及调用示例见 `examples/od_replay.c`。

## 配置和内存

| 字段 | 默认值 | 说明 |
|---|---:|---|
| `observation_capacity` | 30 | 环形位置观测队列点数，2～65536；满后淘汰最旧点 |
| `fit_degree` | 3 | 拟合阶数，1～16，必须小于容量 |
| `nmea_reference_utc_ms` | -1 | 无 RMC 时的 UTC 毫秒参考；-1 表示需要 RMC 日期 |
| `gps_week_rollover` | 2048 | RTCM 1019 的 GPS 周年代偏移，1024 的整数倍，0～8192 |

从 `od_default_config()` 初始化配置后再修改字段。至少 `fit_degree + 1` 个有效观测才可输出，默认第 4 个有效观测开始输出，不必等 30 点队列填满。`od_create(NULL, &ctx)` 同样使用此默认配置，适用于 NMEA 和 RTCM。

在最新观测时间，对队列内全部点做切比雪夫多项式最小二乘拟合，求位置和一阶导数。点间隔可不均匀，时间必须严格递增。默认使用 30 点三次拟合；位置噪声、采样间隔和窗口长度都会影响估速。需要三点二次拟合时，应同时显式设置容量为 3、阶数为 2。

位置队列按容量分配。NMEA 使用固定 1025 字节语句缓冲，RTCM 使用固定 1029 字节组帧缓冲；RTKLIB 解码器在首次非空 RTCM 输入时创建，并另行保留星历与当前观测历元。因此 RTCM 模式的总内存还包括位置队列之外的解码与定位状态。reset 会释放 RTCM 解码器，destroy 释放所有实例资源。

## 数据约定

- 每次输入最多 1 MiB，允许分包和粘包。`length == 0` 为无操作；`out_info` 可以为 NULL。输入缓冲只需在调用期间有效。
- NMEA 带校验和语句在 `*HH` 收齐时处理，无校验和语句需要 CR 或 LF 结束。带校验和时必须校验通过。
- `OD_OK` 表示本次输入已处理或缓存，不保证已有新的位置点。`accepted_observations` 是本次新增位置点数，`rejected_records` 是丢弃记录数，`buffered_observations` 是当前队列点数。
- 坏校验、非法位置、缺日期、重复或倒序位置点会被丢弃并计数，后续有效数据继续处理。等待星历或不支持的消息可能不产生任何点。
- 一个实例只处理一个目标的一种输入来源。切换 RTCM/NMEA 前 reset，或用不同实例。不会自动融合两条来源。
- 同一实例由调用方串行调用。由于 RTKLIB 共享内部状态，多个 RTCM 实例也应串行解码，不承诺并发线程安全。
- `od_get_orbit` 返回最新有效观测时间的状态，不使用主机当前时间、不接受未来查询时间。失败或未就绪时保持输出结构体不变。
- 一次输入产生多个历元时，输出接口只提供最新快照，不另设输出队列。需要全部历元时，应按历元输入并读取。

RTCM 路径：RTCM3 观测和星历 → RTKLIB GPS/北斗单点定位 → 带 UTC 时间的 ECEF 位置 → J2000 拟合状态。沿用原项目单频、广播星历等配置，不新增 RTK、PPP 或原始观测级精密定轨。

NMEA 路径：GGA 提供经纬度、高程，RMC 提供日期。海拔与大地水准面差相加得到椭球高，再转 ECEF；水准面差缺省时按原实现使用 0。只有 RMC 不产生三维位置点。

只有 GGA 时，`nmea_reference_utc_ms` 应与第一条 GGA 相差小于 12 小时；后续可连续跨午夜。间隔达到 12 小时以上，需要新的 RMC，或重建实例并提供新的时间参考。GGA 自身无法分辨跨多天的同一时刻。

RTCM 时间来自 1019/1042 星历和 MSM 历元；默认 GPS 周范围为 2048～3071，其他年代需设置相应偏移。不读取主机日期猜测历史数据。

## 预编译库

Linux 安装包在 `dist/linux-x86_64`，包括 `lib/libod.a`、对外头文件、单文件源码和许可。外部程序同样只需 C 编译器：

```sh
gcc -std=c99 -I/path/to/od/include app.c /path/to/od/lib/libod.a \
    -lm -o app
```

CMake 调用方式：`find_package(od CONFIG REQUIRED)`，然后 `target_link_libraries(app PRIVATE od::orbit_determination)`。独立示例：

```sh
cmake -S examples/prebuilt -B build-consumer \
    -Dod_DIR="$PWD/dist/linux-x86_64/lib/cmake/od"
cmake --build build-consumer
```

交叉编译工具链文件保留在 `cmake/toolchains/riscv64-linux-gnu.cmake`，已去掉 C++ 编译器配置。**按本次要求，没有构建或验证 RISC-V 版本。** 旧的 C++ 版本安装包已归档，不作为本次纯 C 版本交付。

## 验证与来源

基础测试也使用纯 C，覆盖接口、分片、日期、环形缓存、数值拟合及禁止外推。可在配置时增加 `-DOD_RTCM_TEST_FILE=/path/to/capture.rtcm3` 启用真实 RTCM 分片/reset 回放测试。验证详情见 `验证记录.md`。

坐标转换沿用原项目，EOP 2026 数据表已内嵌到实现文件。表外使用端点，TT-UTC 沿用 69.184 秒；其他年代或更高精度任务需要更新模型数据。已有回放回归不等同于独立真值精度评估。

内嵌 RTKLIB 2.4.3 b34 的源代码与 BSD 许可保留原始版权声明。为合并成一个编译单元，文件内私有符号增加了模块前缀，并通过宏保存/恢复隔离各原始源码段的宏定义。原第三方许可和来源说明也保留在 `third_party/rtklib`。

转换前源码和安装包备份保存在 `build-single-c/previous_sources.zip`；旧的多文件实现和旧 RISC-V 安装包已移入 `build-single-c/previous-files`。这些目录属于构建工作区，不参与编译或安装。
