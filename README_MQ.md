# 轨道预报服务 CCSM POSIX MQ 二进制协议

服务不使用 gRPC 或 Protobuf。每次 `mq_send()` 投递一条完整的 CCSM 固定头帧，
接收方使用 `mq_receive()` 返回的实际长度校验帧头中的 `payload_length`。

共享的协议定义和安全编解码实现位于：

- `include/orbit_mq_protocol.h`
- `src/orbit_mq_protocol.cpp`

## 队列

- 生命周期 → 轨道预报：`/csm_main_to_end0`
- 轨道预报 → 生命周期：`/csm_end0_to_main`
- 轨道预报主动心跳/遥测：`/csm_end0_to_main`
- 队列深度：10
- 现场单条 MQ 消息上限：2048 字节；重新创建队列时可以配置为 8192 字节

## 10 字节帧头

多字节数值统一使用小端序。不要把偏移 1、6 强制转换为
`unsigned int *`，U54 上可能产生非对齐访问；应逐字节读写或使用本项目编解码函数。

| 偏移 | 长度 | 类型 | 字段 |
|---:|---:|---|---|
| 0 | 1 | `u8` | `message_type` |
| 1 | 4 | `u32 LE` | `message_seq` |
| 5 | 1 | `u8` | `stream_end` |
| 6 | 4 | `u32 LE` | `payload_length` |
| 10 | N | bytes | `payload` |

总长度必须满足：

```text
mq_receive返回长度 == 10 + payload_length
```

`stream_end`：

```text
0x00 = 后面还有分片
0x01 = 现场 CCSM 兼容结束值（接收端兼容）
0xFF = 单包消息或最后一包
```

轨道预报服务发送时统一使用 `0xFF`；接收时同时接受现场当前使用的
`0x01` 和接口文档规定的 `0xFF`。

参考资料 `mqueue_frame.txt` 的文字把消息长度写成了 `u8`，但
`CCSM_APP_Process.cpp` 实际按 `unsigned int` 读写该字段，因此本实现按
`u32`、10 字节帧头处理。现场读取到的队列属性为 `mq_msgsize=2048`；程序
打开已有队列后会通过 `mq_getattr()` 使用实际上限，因此也兼容配置为 8192
字节的新队列。

## 生命周期 → 轨道预报类型

| 值 | 枚举 | payload |
|---:|---|---|
| 0 | `Unspecified` | 无效/预留 |
| 1 | `ReceiveUplinkData` | `UplinkPacket` |
| 2 | `ReceiveTimeSync` | 忽略；不解析 payload，不改变任何状态 |
| 3 | `ReceiveRtcmData` | 原始 RTCM3 字节 |
| 4 | `GetStatus` | 空 |
| 5 | `Stop` | 空 |
| 6 | `Reset` | 空 |
| 7 | `PredictOrbit` | `OrbitPredictionRequest` |

类型 1、2、3、5、6 为异步消息，轨道服务不发送通用回复。类型 4 返回状态。
类型 7 的预报结果通过类型 4 `ReceiveTelemetryData` 分片输出。

## 轨道预报 → 生命周期类型

| 值 | 枚举 | payload |
|---:|---|---|
| 0 | `Unspecified` | 无效/预留 |
| 1 | `Heartbeat` | `StatusReply`；也用于 `GetStatus` 结果 |
| 2 | `TelemetryData` | 遥测预留 |
| 3 | `TelemetryStatus` | 状态预留 |
| 4 | `ReceiveTelemetryData` | 原始 `OrbitPoint[]`；轨道预报结果 |

## 业务 payload

### UplinkPacket / TelemetryPacket

```text
u32 port
u32 packet_header_length
u8  packet_header[packet_header_length]
u32 data_length
u8  data[data_length]
```

### TimeSyncData

```text
i64 timestamp_s
```

单位为 Unix UTC 秒。当前服务明确忽略该消息，不解析 payload，也不把它用于
RTCM 解码或轨道预报。RTCM 时间完全由 1019/1042 星历和 MSM 历元建立。

### RTCMData

payload 全部是原始 RTCM3 字节，不再增加内部长度。每包最大负载根据
`mq_getattr()` 得到的实际 `mq_msgsize - 10` 计算：现场 2048 字节队列对应
2038 字节；8192 字节队列对应 8182 字节。`message_seq` 从 0 递增，中间包
`stream_end=0x00`，最后一包发送 `0xFF`，并兼容接收现场的 `0x01`。服务端
RTCM 帧缓存支持 RTCM3 帧跨 MQ 分片。

### OrbitPredictionRequest

```text
i64 start_time_s
u32 duration_s
```

共 12 字节，时间单位为 Unix UTC 秒。

> **临时现场测试逻辑：** 服务暂时忽略 `start_time_s`，强制使用最近一次 RTCM
> ECEF 解算时间作为预报起点。完成 RTCM/MQ 联调后应删除该覆盖逻辑。

### StatusReply / Heartbeat

```text
i64 timestamp_ms
u8  service_status       # 0未知、1空闲、2运行、3异常
u32 message_length
u8  message[message_length]  # UTF-8
```

心跳服务启动后立即发送一次，之后默认每 1000 ms 发送一次。心跳是尽力投递，
不等待回复，队列不存在或已满不会使轨道服务退出。

### OrbitData

类型 4 的 MQ payload 直接由连续的 `OrbitPoint[]` 构成；没有
`TelemetryPacket` 外壳，也没有 `point_count`。接收端通过
`payload_length / 56` 得到本帧点数，payload 长度必须能被 56 整除。

每个 `OrbitPoint` 固定 56 字节：

```text
i64 timestamp_ms
f64 x
f64 y
f64 z
f64 vx
f64 vy
f64 vz
```

轨道结果帧的 `message_seq` 从 0 连续递增。数据帧使用
`stream_end=0x00`；最后额外发送一条 payload 为空、`stream_end=0xFF` 的
终止帧。接收方应检查序号连续性，终止帧表示本次输出结束。

## 配置

```text
ORBIT_MQ_REQUEST_QUEUE=/csm_main_to_end0
ORBIT_MQ_RESPONSE_QUEUE=/csm_end0_to_main
ORBIT_MQ_MAX_MESSAGES=10
ORBIT_MQ_MESSAGE_SIZE=2048
ORBIT_MQ_CREATE=0
ORBIT_LIFECYCLE_REQUEST_QUEUE=/csm_end0_to_main
ORBIT_HEARTBEAT_ENABLED=1
ORBIT_HEARTBEAT_INTERVAL_MS=1000
ORBIT_CONSOLE_LOG=1
ORBIT_OBSERVATION_CAPACITY=600
ORBIT_FIT_DEGREE=10
ORBIT_STREAM_BATCH_SIZE=100
ORBIT_RTCM_CACHE_BYTES=4194304
```

## Ubuntu 24.04 chroot / GCC 13 静态构建

RTKLIB 2.4.3 源码直接编译进程序。MQ 传输不再依赖 Protobuf、gRPC、
OpenSSL、Abseil、c-ares 或 RE2。

```bash
chroot /root/.cache/orbit-noble-amd64 /bin/bash -lc '
  cd /mnt/e/JC/OrbitDetermination/PredictOrbit
  bash scripts/build_riscv64_static.sh
'
```

产物：

```text
build-riscv64-static/orbit_prediction_server_static
build-riscv64-static/orbit_mq_cli_static
```

脚本校验完全静态链接和 U54 `rv64gc/lp64d` ISA，并拒绝 RVV/B 扩展。

## Windows 可视化测试

```powershell
powershell -ExecutionPolicy Bypass -File tools\run_windows_mq_tester.ps1
```

Windows 工具通过 SSH 上传静态服务端和 CLI。RTCM 抓取会先清洗 CRC 有效帧，
自动授时，然后使用上述固定二进制帧发送；预报输出下载为 CSV 并显示轨迹。

## 测试

`orbit_mq_protocol_test` 校验逐字节帧布局、长度校验及业务 payload 往返；
`orbit_prediction_mq_smoke_test` 覆盖 7 个入站类型、状态返回和预报分片。

### 192.168.104.100 现场 RTCM 验证（2026-08-27）

- `/csm_main_to_end0`：`mq_maxmsg=10`，`mq_msgsize=2048`
- 实际帧：`type=3`、`sequence=0`、`stream_end=0x01`
- 抽样消息：总长 732 字节，payload 722 字节，长度字段完全匹配
- 抽样 payload 含 6 个 CRC24Q 正确的 RTCM3 帧，无损坏帧或前导杂字节
- 临时启动本版本服务后共接收并解析 60 个 RTCM3 帧，队列积压降为 0

现场这一小批数据未包含足够的连续星历和观测信息，因而尚未产生 ECEF
定位结果；需要上游继续发送完整 RTCM 数据段后再进行定位和 3600 秒预测测试。
抽样原始固定帧保存在 `out/live_rtcm_frame_20260814.bin`。
