cmake_minimum_required(VERSION 3.16)

project(calgnss_orbit_mini C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CALGNSS_ORBIT_MINI_SOURCES
    src/calgnss_orbit_mini.c
)

add_library(calgnss_orbit_mini SHARED
    ${CALGNSS_ORBIT_MINI_SOURCES}
)

add_library(calgnss_orbit_mini_static STATIC
    ${CALGNSS_ORBIT_MINI_SOURCES}
)

set_target_properties(calgnss_orbit_mini PROPERTIES
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set_target_properties(calgnss_orbit_mini_static PROPERTIES
    OUTPUT_NAME calgnss_orbit_mini_static
)

target_include_directories(calgnss_orbit_mini PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_include_directories(calgnss_orbit_mini_static PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (UNIX)
    target_link_libraries(calgnss_orbit_mini PUBLIC m)
    target_link_libraries(calgnss_orbit_mini_static PUBLIC m)
endif()
