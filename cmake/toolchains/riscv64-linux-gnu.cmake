set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR riscv64)

set(CMAKE_C_COMPILER riscv64-linux-gnu-gcc)
set(CMAKE_CXX_COMPILER riscv64-linux-gnu-g++)
set(CMAKE_AR riscv64-linux-gnu-ar)
set(CMAKE_RANLIB riscv64-linux-gnu-ranlib)
set(CMAKE_STRIP riscv64-linux-gnu-strip)

# A conservative Linux baseline that works on considerably more boards than
# the Ubuntu compiler's distribution-specific default. Override these cache
# variables for a known CPU (for example rv64gcv) when appropriate.
set(ORBIT_RISCV_MARCH "rv64gc" CACHE STRING "RISC-V ISA used for target objects")
set(ORBIT_RISCV_MABI "lp64d" CACHE STRING "RISC-V GNU/Linux ABI")
string(APPEND CMAKE_C_FLAGS_INIT " -march=${ORBIT_RISCV_MARCH} -mabi=${ORBIT_RISCV_MABI}")
string(APPEND CMAKE_CXX_FLAGS_INIT " -march=${ORBIT_RISCV_MARCH} -mabi=${ORBIT_RISCV_MABI}")

# Ubuntu/Debian multiarch target packages install below /usr while the cross
# compiler runtime and headers also live below /usr/riscv64-linux-gnu.
set(CMAKE_FIND_ROOT_PATH
    /usr/riscv64-linux-gnu
    /usr
)
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Enables CTest for binaries which only depend on the cross sysroot. Tests that
# need additional target shared libraries may also require QEMU_LD_PREFIX=/.
set(CMAKE_CROSSCOMPILING_EMULATOR /usr/bin/qemu-riscv64;-L;/usr/riscv64-linux-gnu)
