cmake_minimum_required(VERSION 3.16)

project(CalGNSSC C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(calgnss_orbit_mini STATIC
    src/calgnss_orbit_mini.c
)

target_include_directories(calgnss_orbit_mini PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (UNIX)
    target_link_libraries(calgnss_orbit_mini PUBLIC m)
endif()

add_executable(CalGNSSC
    examples/export_j2000.c
)

target_link_libraries(CalGNSSC PRIVATE calgnss_orbit_mini)
