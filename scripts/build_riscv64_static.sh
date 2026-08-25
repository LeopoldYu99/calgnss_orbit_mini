#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${project_dir}/build-riscv64-static"
prefix="${ORBIT_STATIC_PREFIX:-/root/.cache/orbit-riscv64-static/prefix}"

if [[ ! -f "${prefix}/lib/libgrpc++.a" || ! -f "${prefix}/lib/libprotobuf.a" ]]; then
    echo "Static RISC-V dependencies are missing from ${prefix}." >&2
    echo "Run scripts/build_riscv64_static_deps.sh first." >&2
    exit 1
fi

export PKG_CONFIG_LIBDIR="${prefix}/lib/pkgconfig"
export PKG_CONFIG_SYSROOT_DIR="/"

cmake --fresh -S "${project_dir}" -B "${build_dir}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE="${project_dir}/cmake/toolchains/riscv64-linux-gnu-static.cmake" \
    -DORBIT_RISCV_STATIC_PREFIX="${prefix}" \
    -DCMAKE_PREFIX_PATH="${prefix}" \
    -DProtobuf_DIR="${prefix}/lib/cmake/protobuf" \
    -DgRPC_DIR="${prefix}/lib/cmake/grpc" \
    -Dabsl_DIR="${prefix}/lib/cmake/absl" \
    -Dre2_DIR="${prefix}/lib/cmake/re2" \
    -Dc-ares_DIR="${prefix}/lib/cmake/c-ares" \
    -DOPENSSL_USE_STATIC_LIBS=TRUE \
    -DOPENSSL_ROOT_DIR="${prefix}" \
    -DZLIB_ROOT="${prefix}" \
    -DCALGNSS_BUILD_GRPC=ON \
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
