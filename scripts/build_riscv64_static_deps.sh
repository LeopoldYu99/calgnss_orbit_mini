#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_root="${ORBIT_STATIC_SOURCE_ROOT:-/root/.cache/orbit-riscv64-static/ubuntu-sources}"
work_root="${ORBIT_STATIC_WORK_ROOT:-/root/.cache/orbit-riscv64-static/build}"
prefix="${ORBIT_STATIC_PREFIX:-/root/.cache/orbit-riscv64-static/prefix}"
toolchain="${project_dir}/cmake/toolchains/riscv64-linux-gnu-static.cmake"
parallel="${ORBIT_BUILD_PARALLEL:-$(nproc)}"

find_source() {
    local pattern="$1"
    local matches=("${source_root}"/${pattern})
    if [[ ! -d "${matches[0]}" ]]; then
        echo "Missing source ${source_root}/${pattern}" >&2
        echo "Enable Ubuntu deb-src entries and download the required source packages first." >&2
        return 1
    fi
    printf '%s\n' "${matches[0]}"
}

cmake_static() {
    local name="$1"
    local source_dir="$2"
    shift 2
    cmake --fresh -S "${source_dir}" -B "${work_root}/${name}" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_TOOLCHAIN_FILE="${toolchain}" \
        -DORBIT_RISCV_STATIC_PREFIX="${prefix}" \
        -DCMAKE_INSTALL_PREFIX="${prefix}" \
        -DCMAKE_PREFIX_PATH="${prefix}" \
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
        -DBUILD_SHARED_LIBS=OFF \
        "$@"
    cmake --build "${work_root}/${name}" --parallel "${parallel}"
    cmake --install "${work_root}/${name}"
}

mkdir -p "${work_root}" "${prefix}"

# Ubuntu 26.04 and newer may ship RISC-V static runtime objects built for a
# distribution baseline that includes V/B extensions. Project-level
# -march=rv64gc cannot downgrade instructions already present in libc.a.
# Reject that toolchain before rebuilding all dependencies.
probe_file="${work_root}/toolchain-isa-probe"
printf 'int main(void) { return 0; }\n' | \
    riscv64-linux-gnu-gcc -x c - -O2 -static -march=rv64gc -mabi=lp64d \
        -o "${probe_file}"
probe_arch="$(riscv64-linux-gnu-readelf -A "${probe_file}" | \
    sed -n 's/.*Tag_RISCV_arch: "\([^"]*\)"/\1/p')"
if [[ "${probe_arch}" == *"_v"* || "${probe_arch}" == *"_zv"* ||
      "${probe_arch}" == *"_b"* || "${probe_arch}" == *"_zba"* ||
      "${probe_arch}" == *"_zbb"* || "${probe_arch}" == *"_zbs"* ]]; then
    echo "Static runtime ISA is not compatible with an rv64gc/U54 target:" >&2
    echo "  ${probe_arch}" >&2
    echo "Use an Ubuntu 24.04 GCC 13 cross environment or the board SDK." >&2
    exit 1
fi
echo "Static runtime ISA check passed: ${probe_arch}"

zlib_source="$(find_source 'zlib-*')"
openssl_source="$(find_source 'openssl-*')"
absl_source="$(find_source 'abseil-*')"
protobuf_source="$(find_source 'protobuf-*')"
cares_source="$(find_source 'c-ares-*')"
re2_source="$(find_source 're2-*')"
grpc_source="$(find_source 'grpc-*')"

if [[ ! -f "${zlib_source}/zconf.h" && -f "${zlib_source}/zconf.h.included" ]]; then
    cmake -E rename "${zlib_source}/zconf.h.included" "${zlib_source}/zconf.h"
fi
cmake -E make_directory "${work_root}/zlib-autotools"
pushd "${work_root}/zlib-autotools" >/dev/null
CHOST=riscv64-linux-gnu \
CC=riscv64-linux-gnu-gcc \
AR=riscv64-linux-gnu-ar \
RANLIB=riscv64-linux-gnu-ranlib \
CFLAGS="-O3 -march=rv64gc -mabi=lp64d" \
"${zlib_source}/configure" --static --prefix="${prefix}"
make -j"${parallel}"
make install
popd >/dev/null

cmake -E make_directory "${work_root}/openssl"
pushd "${work_root}/openssl" >/dev/null
CC=riscv64-linux-gnu-gcc \
CXX=riscv64-linux-gnu-g++ \
AR=riscv64-linux-gnu-ar \
RANLIB=riscv64-linux-gnu-ranlib \
"${openssl_source}/Configure" linux64-riscv64 \
    no-shared no-tests no-apps no-docs no-module \
    -march=rv64gc -mabi=lp64d \
    --prefix="${prefix}" --libdir=lib
make -j"${parallel}"
make install_sw
popd >/dev/null

cmake_static abseil "${absl_source}" \
    -DABSL_BUILD_TESTING=OFF \
    -DABSL_ENABLE_INSTALL=ON \
    -DCMAKE_CXX_STANDARD=17

cmake_static protobuf "${protobuf_source}" \
    -Dprotobuf_BUILD_TESTS=OFF \
    -Dprotobuf_BUILD_CONFORMANCE=OFF \
    -Dprotobuf_BUILD_EXAMPLES=OFF \
    -Dprotobuf_BUILD_PROTOC_BINARIES=OFF \
    -Dprotobuf_BUILD_LIBPROTOC=ON \
    -Dprotobuf_BUILD_SHARED_LIBS=OFF \
    -Dprotobuf_WITH_ZLIB=ON \
    -DZLIB_ROOT="${prefix}"

cmake_static cares "${cares_source}" \
    -DCARES_SHARED=OFF \
    -DCARES_STATIC=ON \
    -DCARES_BUILD_TESTS=OFF \
    -DCARES_BUILD_TOOLS=OFF

cmake_static re2 "${re2_source}" \
    -DRE2_BUILD_TESTING=OFF \
    -DRE2_TEST=OFF \
    -DRE2_BENCHMARK=OFF

cmake_static grpc "${grpc_source}" \
    -DgRPC_INSTALL=ON \
    -DgRPC_BUILD_TESTS=OFF \
    -DgRPC_BUILD_CODEGEN=OFF \
    -DgRPC_BUILD_GRPC_CPP_PLUGIN=OFF \
    -DgRPC_BUILD_GRPC_CSHARP_PLUGIN=OFF \
    -DgRPC_BUILD_GRPC_NODE_PLUGIN=OFF \
    -DgRPC_BUILD_GRPC_OBJECTIVE_C_PLUGIN=OFF \
    -DgRPC_BUILD_GRPC_PHP_PLUGIN=OFF \
    -DgRPC_BUILD_GRPC_PYTHON_PLUGIN=OFF \
    -DgRPC_BUILD_GRPC_RUBY_PLUGIN=OFF \
    -DgRPC_ABSL_PROVIDER=package \
    -DgRPC_CARES_PROVIDER=package \
    -DgRPC_PROTOBUF_PROVIDER=package \
    -DgRPC_RE2_PROVIDER=package \
    -DgRPC_SSL_PROVIDER=package \
    -DgRPC_ZLIB_PROVIDER=package \
    -DOPENSSL_USE_STATIC_LIBS=TRUE \
    -DOPENSSL_ROOT_DIR="${prefix}" \
    -DProtobuf_PROTOC_EXECUTABLE=/usr/bin/protoc

echo "Static RISC-V dependencies installed in ${prefix}"
