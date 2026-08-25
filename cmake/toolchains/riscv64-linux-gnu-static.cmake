include(${CMAKE_CURRENT_LIST_DIR}/riscv64-linux-gnu.cmake)

set(ORBIT_RISCV_STATIC_PREFIX "" CACHE PATH
    "Prefix containing statically built RISC-V dependencies")
if (ORBIT_RISCV_STATIC_PREFIX)
    list(PREPEND CMAKE_FIND_ROOT_PATH "${ORBIT_RISCV_STATIC_PREFIX}")
endif()

set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
string(APPEND CMAKE_EXE_LINKER_FLAGS_INIT " -static")

