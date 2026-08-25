# OrbitPredictionService gRPC 微服务

## 服务信息

- Proto package：`orbit_prediction`
- Service：`orbit_prediction.OrbitPredictionService`
- 生产地址：`192.168.104.100:50051`
- 传输：当前使用内网明文 gRPC；跨不可信网络部署前应启用 TLS 和鉴权。

服务监听地址由 `ORBIT_GRPC_ADDRESS` 配置。生产 systemd 配置使用
`192.168.104.100:50051`。本机、WSL2 或容器测试应使用
`ORBIT_GRPC_ADDRESS=0.0.0.0:50051`，再从客户端访问 `localhost:50051`。

## RPC 行为

### ReceiveTimeSync

`timestamp_ms` 是 Unix epoch UTC 毫秒。时间同步为只有 UTC 时分秒的 GGA 数据提供日期，
并自动处理最接近同步时间的前后午夜。RMC 中的有效日期优先于时间同步日期。

### ReceiveUplinkData

`data` 是 ASCII/UTF-8 NMEA 文本，可以一次包含多行。当前支持：

- RMC：读取 UTC 日期，状态必须为 `A` 才会采用。
- GGA：读取 UTC 时间、经纬度、定位质量、MSL 高和 geoid separation，并转换为
  WGS84 ECEF 米坐标后写入算法缓存。

GGA 的 fix quality 必须大于 0；观测必须严格按时间递增。如果句子包含 `*HH` 校验和，
服务会验证校验和。`packet_header` 和 `port` 保留给上行传输层，当前不参与轨道计算。
解析或校验失败时 RPC 本身正常返回，但 `CommonReply.success=false`，服务日志会记录原因。

### ReceiveRTCMData

`data` 必须是 RTCM3 原始二进制字节，不是十六进制文本。接口支持：

- 一个 RPC 中携带一个或多个 RTCM3 帧。
- 一个 RTCM3 帧拆分到多次 RPC（半包）。
- 帧前噪声和粘包重同步。
- RTCM3 长度字段和 CRC24Q 校验。
- 解析消息类型并返回在 `CommonReply.message` 中。
- 默认缓存最近 4 MiB 的完整帧，可通过 `ORBIT_RTCM_CACHE_BYTES` 调整。

RTCM MSM（例如附件中的 1077、1097、1117、1127）携带原始伪距/载波相位，
1019/1042 携带星历；它们不是接收机 ECEF 位置。当前接口负责接收、校验和缓存，
不会把 RTCM 帧直接写入轨道拟合缓存。要让 RTCM 独立驱动轨道预报，还需要接入
RTKLIB 或等价 GNSS PVT 解算器，先从观测和星历求出带 UTC 的 ECEF 位置。

### PredictOrbit

至少通过 NMEA/PVT 路径得到两组有效 ECEF 观测后才能调用。请求参数为：

- `start_time_s`：Unix UTC 秒，输出包含这个起始点。
- `duration_s`：持续时间，范围 1～3600 秒。
- `step_s`：步长，范围 0.001～`duration_s` 秒，支持小数。

输出点不超过 100000 个，每个 `OrbitData` 默认包含 10 个 `OrbitPoint`，批大小通过
`ORBIT_STREAM_BATCH_SIZE` 修改。输出坐标系为 J2000，位置单位为米，速度单位为米/秒。
起始时间不得早于当前观测缓存，结束时间不得晚于最新观测后 3600 秒。

### Stop 和 Reset

`Stop` 终止活动预测流，并阻止新的 `PredictOrbit` 调用，但不会关闭 gRPC 进程。
`Reset` 清空观测、RTCM 缓存、时间同步和 Stop 状态，之后可以重新接收数据和预测。

所有普通 RPC 都通过 `CommonReply.success` 表示结果，并在 `CommonReply.message` 中返回
成功摘要或失败原因。流式 `PredictOrbit` 的请求错误通过 gRPC status 返回。

## WSL2 构建

管理员 PowerShell：

```powershell
wsl --install -d Ubuntu
```

重启 Windows 后进入 Ubuntu。建议将仓库放在 WSL 的 Linux 文件系统以获得更好的构建性能；
直接使用当前 Windows 工作区时路径为：

```bash
cd /mnt/e/JC/OrbitDetermination/PredictOrbit
bash scripts/setup_wsl_ubuntu.sh
```

不要复用 Windows 的 `build` 目录；脚本使用独立的 `build-wsl`。

## RISC-V 64 交叉编译

已经提供 `riscv64-linux-gnu` CMake 工具链。构建机需要宿主架构的 `protoc`、
`grpc_cpp_plugin`，以及目标架构的编译器、gRPC、Protobuf 和 OpenSSL 等开发包。
当前 WSL Ubuntu 已使用 multiarch 安装这些依赖。执行：

```bash
cd /mnt/e/JC/OrbitDetermination/PredictOrbit
bash scripts/build_riscv64.sh
```

输出文件为 `build-riscv64/orbit_prediction_server`，默认使用较通用的
`rv64gc` 指令集和 `lp64d` ABI。已知目标 CPU 支持其他扩展时可以覆盖，例如：

```bash
cmake --fresh -S . -B build-riscv64 -G Ninja \
  -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/riscv64-linux-gnu.cmake \
  -DORBIT_RISCV_MARCH=rv64gcv -DORBIT_RISCV_MABI=lp64d
cmake --build build-riscv64 --parallel
```

在非 RISC-V 构建机上可使用 QEMU 验证：

```bash
ORBIT_GRPC_ADDRESS=127.0.0.1:50052 \
qemu-riscv64 -L /usr/riscv64-linux-gnu \
  -E LD_LIBRARY_PATH=/usr/lib/riscv64-linux-gnu \
  ./build-riscv64/orbit_prediction_server
```

另一个终端生成 Python stub 并测试全部 RPC：

```bash
mkdir -p build-riscv64/generated-python
protoc -I proto -I /usr/include \
  --python_out=build-riscv64/generated-python \
  --grpc_out=build-riscv64/generated-python \
  --plugin=protoc-gen-grpc=/usr/bin/grpc_python_plugin \
  proto/orbit_prediction.proto
PYTHONPATH=build-riscv64/generated-python \
  python3 tests/grpc_smoke_test.py --address 127.0.0.1:50052
```

### RISC-V 全静态版本

Icicle Kit 的 SiFive U54 只支持 `rv64imafdc`（`rv64gc`），不支持 RISC-V
V/B 扩展。请使用 Ubuntu 24.04 的 GCC 13 交叉工具链或板卡配套 SDK 构建。
Ubuntu 26.04 的 RISC-V 静态运行库可能已经包含 V/B 指令；即使项目源码指定
`-march=rv64gc`，链接后的程序仍会在 U54 上触发 `illegal instruction`。构建脚本
会检查静态运行库和最终 ELF 的 ISA 属性并拒绝不兼容的产物。

全静态构建会从源码生成 gRPC、Protobuf、OpenSSL、Abseil、RE2、c-ares 和 zlib 的
RISC-V `.a` 库。Ubuntu 源码包准备完成后执行：

```bash
bash scripts/build_riscv64_static_deps.sh
bash scripts/build_riscv64_static.sh
```

部署文件为 `build-riscv64-static/orbit_prediction_server_static`。它是已裁剪的
`rv64gc/lp64d` 全静态 ELF，不需要在目标机安装 gRPC、Protobuf 或 C++ 运行库：

```bash
chmod +x orbit_prediction_server_static
ORBIT_GRPC_ADDRESS=0.0.0.0:50051 ./orbit_prediction_server_static
```

构建目录同时保留未裁剪的 `orbit_prediction_server`，供符号分析和故障排查。

Windows 原生构建也可使用 Visual Studio 自带的 vcpkg manifest：

```powershell
cmake -S . -B build-grpc `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake"
cmake --build build-grpc --config Release
```

## 运行

WSL2 测试：

```bash
ORBIT_GRPC_ADDRESS=0.0.0.0:50051 ./build-wsl/orbit_prediction_server
```

如果目标 Linux 主机确实配置了 `192.168.104.100`：

```bash
ORBIT_GRPC_ADDRESS=192.168.104.100:50051 ./build-wsl/orbit_prediction_server
```

检查地址是否存在：

```bash
ip address show | grep 192.168.104.100
```

## systemd 部署

```bash
sudo useradd --system --no-create-home --shell /usr/sbin/nologin orbit
sudo install -d -o root -g root /opt/orbit-prediction/bin
sudo install -m 0755 build-wsl/orbit_prediction_server /opt/orbit-prediction/bin/
sudo install -m 0644 deploy/orbit-prediction.env /etc/orbit-prediction.env
sudo install -m 0644 deploy/orbit-prediction.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now orbit-prediction
sudo systemctl status orbit-prediction
```

## Docker

```bash
docker build -t orbit-prediction .
docker run --rm -p 50051:50051 orbit-prediction
```

此时客户端连接 `localhost:50051`。生产主机需要暴露指定地址时，可使用：

```bash
docker run --rm -p 192.168.104.100:50051:50051 orbit-prediction
```
