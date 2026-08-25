#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${project_dir}/build-riscv64"

# Keep pkg-config from leaking host (x86_64) libraries into the target link.
export PKG_CONFIG_LIBDIR="/usr/lib/riscv64-linux-gnu/pkgconfig:/usr/share/pkgconfig"
export PKG_CONFIG_SYSROOT_DIR="/"

cmake --fresh -S "${project_dir}" -B "${build_dir}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE="${project_dir}/cmake/toolchains/riscv64-linux-gnu.cmake" \
    -DCMAKE_PREFIX_PATH="/usr/lib/riscv64-linux-gnu/cmake" \
    -DOPENSSL_ROOT_DIR="/usr" \
    -DOPENSSL_CRYPTO_LIBRARY="/usr/lib/riscv64-linux-gnu/libcrypto.so" \
    -DOPENSSL_SSL_LIBRARY="/usr/lib/riscv64-linux-gnu/libssl.so" \
    -DCALGNSS_BUILD_GRPC=ON \
    -DCALGNSS_BUILD_TESTS=ON

cmake --build "${build_dir}" --parallel
ctest --test-dir "${build_dir}" --output-on-failure

file "${build_dir}/orbit_prediction_server"
