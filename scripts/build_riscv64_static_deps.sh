#!/usr/bin/env bash
set -euo pipefail

# The fixed-binary MQ transport has no separately built third-party dependency.
# RTKLIB is compiled directly into the application. This compatibility script
# only verifies the Ubuntu 24.04 GCC 13 static sysroot.
work_root="${ORBIT_STATIC_WORK_ROOT:-/root/.cache/orbit-riscv64-static/build}"
mkdir -p "${work_root}"

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
    exit 1
fi
echo "Static runtime ISA check passed: ${probe_arch}"
echo "No Protobuf/gRPC transport dependencies need to be built."
