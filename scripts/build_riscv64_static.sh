#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${project_dir}/build-riscv64-static"
source_root="${ORBIT_STATIC_SOURCE_ROOT:-/root/.cache/orbit-riscv64-static/ubuntu-sources}"
rtklib_source="${ORBIT_RTKLIB_SOURCE_DIR:-${source_root}/rtklib}"

if [[ ! -f "${rtklib_source}/src/rtklib.h" ]]; then
    echo "RTKLIB 2.4.3 source is missing from ${rtklib_source}." >&2
    echo "Set ORBIT_RTKLIB_SOURCE_DIR to the Ubuntu 24.04 RTKLIB source tree." >&2
    exit 1
fi

cmake --fresh -S "${project_dir}" -B "${build_dir}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE="${project_dir}/cmake/toolchains/riscv64-linux-gnu-static.cmake" \
    -DORBIT_RTKLIB_SOURCE_DIR="${rtklib_source}" \
    -DCALGNSS_BUILD_MQ=ON \
    -DCALGNSS_BUILD_TESTS=ON

cmake --build "${build_dir}" --parallel
ctest --test-dir "${build_dir}" --output-on-failure

if command -v file >/dev/null 2>&1; then
    file "${build_dir}/orbit_prediction_server"
fi
if riscv64-linux-gnu-readelf -d "${build_dir}/orbit_prediction_server" | grep -q NEEDED; then
    echo "Unexpected dynamic dependency in static server" >&2
    exit 1
fi

# Keep the unstripped executable for diagnostics and create the smaller file
# intended to be copied directly to the target board.
cmake -E copy "${build_dir}/orbit_prediction_server" \
    "${build_dir}/orbit_prediction_server_static"
riscv64-linux-gnu-strip --strip-unneeded \
    "${build_dir}/orbit_prediction_server_static"
cmake -E copy "${build_dir}/orbit_mq_cli" \
    "${build_dir}/orbit_mq_cli_static"
riscv64-linux-gnu-strip --strip-unneeded \
    "${build_dir}/orbit_mq_cli_static"
static_arch="$(riscv64-linux-gnu-readelf -A \
    "${build_dir}/orbit_prediction_server_static" | \
    sed -n 's/.*Tag_RISCV_arch: "\([^"]*\)"/\1/p')"
if [[ "${static_arch}" == *"_v"* || "${static_arch}" == *"_zv"* ||
      "${static_arch}" == *"_b"* || "${static_arch}" == *"_zba"* ||
      "${static_arch}" == *"_zbb"* || "${static_arch}" == *"_zbs"* ]]; then
    echo "Final executable requires instructions outside the rv64gc/U54 baseline:" >&2
    echo "  ${static_arch}" >&2
    exit 1
fi
echo "Final executable ISA: ${static_arch}"
if command -v file >/dev/null 2>&1; then
    file "${build_dir}/orbit_prediction_server_static"
fi
sha256sum "${build_dir}/orbit_prediction_server_static"
if riscv64-linux-gnu-readelf -d "${build_dir}/orbit_mq_cli_static" | grep -q NEEDED; then
    echo "Unexpected dynamic dependency in static MQ CLI" >&2
    exit 1
fi
cli_arch="$(riscv64-linux-gnu-readelf -A \
    "${build_dir}/orbit_mq_cli_static" | \
    sed -n 's/.*Tag_RISCV_arch: "\([^"]*\)"/\1/p')"
if [[ "${cli_arch}" == *"_v"* || "${cli_arch}" == *"_zv"* ||
      "${cli_arch}" == *"_b"* || "${cli_arch}" == *"_zba"* ||
      "${cli_arch}" == *"_zbb"* || "${cli_arch}" == *"_zbs"* ]]; then
    echo "MQ CLI requires instructions outside the rv64gc/U54 baseline:" >&2
    echo "  ${cli_arch}" >&2
    exit 1
fi
echo "MQ CLI ISA: ${cli_arch}"
sha256sum "${build_dir}/orbit_mq_cli_static"
