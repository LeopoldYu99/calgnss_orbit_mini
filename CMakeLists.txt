cmake_minimum_required(VERSION 3.16)

project(calgnss_orbit_mini LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CALGNSS_BUILD_MQ "Build the OrbitPredictionService POSIX MQ service" ON)
option(CALGNSS_BUILD_TESTS "Build project tests" ON)
set(ORBIT_RTKLIB_SOURCE_DIR "" CACHE PATH "RTKLIB 2.4.3 source directory")

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
    POSITION_INDEPENDENT_CODE ON
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

find_package(Threads REQUIRED)

add_library(orbit_prediction_engine STATIC
    src/orbit_prediction_engine.cpp
)
target_include_directories(orbit_prediction_engine PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(orbit_prediction_engine PUBLIC calgnss_orbit_mini_static Threads::Threads)
if (MSVC)
    target_compile_options(orbit_prediction_engine PRIVATE /W4 /permissive-)
else()
    target_compile_options(orbit_prediction_engine PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (CALGNSS_BUILD_TESTS)
    enable_testing()
    add_executable(orbit_prediction_engine_test tests/orbit_prediction_engine_test.cpp)
    target_link_libraries(orbit_prediction_engine_test PRIVATE orbit_prediction_engine)
    add_test(NAME orbit_prediction_engine_test COMMAND orbit_prediction_engine_test)
endif()

if (CALGNSS_BUILD_MQ)
    if (NOT UNIX)
        message(FATAL_ERROR "The POSIX MQ service requires a Unix-like target")
    endif()
    if (NOT ORBIT_RTKLIB_SOURCE_DIR OR
        NOT EXISTS "${ORBIT_RTKLIB_SOURCE_DIR}/src/rtklib.h")
        message(FATAL_ERROR
            "ORBIT_RTKLIB_SOURCE_DIR must point to an RTKLIB 2.4.3 source tree")
    endif()
    add_library(orbit_rtklib STATIC
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/rtkcmn.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/rtcm.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/rtcm2.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/rtcm3.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/rtcm3e.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/ephemeris.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/ionex.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/pntpos.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/preceph.c
        ${ORBIT_RTKLIB_SOURCE_DIR}/src/sbas.c
    )
    target_include_directories(orbit_rtklib PUBLIC ${ORBIT_RTKLIB_SOURCE_DIR}/src)
    target_compile_definitions(orbit_rtklib PUBLIC
        NO_STRICT ENAGLO ENAQZS ENAGAL ENACMP ENAIRN NFREQ=5 NEXOBS=3)
    target_link_libraries(orbit_rtklib PUBLIC m Threads::Threads)
    if (NOT MSVC)
        target_compile_options(orbit_rtklib PRIVATE
            -Wno-unused-const-variable -Wno-discarded-qualifiers)
    endif()

    add_library(orbit_rtcm_positioning STATIC service/rtcm_position_solver.cpp)
    target_include_directories(orbit_rtcm_positioning PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/service)
    target_link_libraries(orbit_rtcm_positioning PUBLIC orbit_rtklib)
    target_compile_definitions(orbit_prediction_engine PUBLIC
        ORBIT_ENABLE_RTCM_POSITIONING=1)
    target_include_directories(orbit_prediction_engine PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/service)
    target_link_libraries(orbit_prediction_engine PUBLIC orbit_rtcm_positioning)

    if (CALGNSS_BUILD_TESTS)
        add_executable(rtcm_position_solver_test tests/rtcm_position_solver_test.cpp)
        target_link_libraries(rtcm_position_solver_test PRIVATE orbit_rtcm_positioning)
        add_test(NAME rtcm_position_solver_test COMMAND rtcm_position_solver_test)
    endif()

    add_library(orbit_mq_protocol STATIC src/orbit_mq_protocol.cpp)
    target_include_directories(orbit_mq_protocol PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include)
    if (MSVC)
        target_compile_options(orbit_mq_protocol PRIVATE /W4 /permissive-)
    else()
        target_compile_options(orbit_mq_protocol PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_library(orbit_prediction_mq_service STATIC
        service/heartbeat_mq_publisher.cpp
        service/orbit_prediction_service.cpp
        service/posix_mq_server.cpp
        service/telemetry_mq_publisher.cpp
    )
    target_include_directories(orbit_prediction_mq_service PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/service
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(orbit_prediction_mq_service PUBLIC
        orbit_prediction_engine
        orbit_mq_protocol
        Threads::Threads
        rt
    )
    if (MSVC)
        target_compile_options(orbit_prediction_mq_service PRIVATE /W4 /permissive-)
    else()
        target_compile_options(orbit_prediction_mq_service PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_executable(orbit_prediction_server service/main.cpp)
    target_link_libraries(orbit_prediction_server PRIVATE orbit_prediction_mq_service)
    if (MSVC)
        target_compile_options(orbit_prediction_server PRIVATE /W4 /permissive-)
    else()
        target_compile_options(orbit_prediction_server PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_executable(orbit_mq_cli tools/orbit_mq_cli.cpp)
    target_link_libraries(orbit_mq_cli PRIVATE
        orbit_mq_protocol
        Threads::Threads
        rt
    )
    if (MSVC)
        target_compile_options(orbit_mq_cli PRIVATE /W4 /permissive-)
    else()
        target_compile_options(orbit_mq_cli PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    add_executable(rtcm_position_debug tools/rtcm_position_debug.cpp)
    target_link_libraries(rtcm_position_debug PRIVATE orbit_rtcm_positioning)
    if (MSVC)
        target_compile_options(rtcm_position_debug PRIVATE /W4 /permissive-)
    else()
        target_compile_options(rtcm_position_debug PRIVATE -Wall -Wextra -Wpedantic)
    endif()

    if (CALGNSS_BUILD_TESTS)
        add_executable(orbit_prediction_mq_smoke_test tests/mq_smoke_test.cpp)
        target_link_libraries(orbit_prediction_mq_smoke_test PRIVATE
            orbit_mq_protocol
            Threads::Threads
            rt
        )
        if (MSVC)
            target_compile_options(orbit_prediction_mq_smoke_test PRIVATE /W4 /permissive-)
        else()
            target_compile_options(orbit_prediction_mq_smoke_test PRIVATE
                -Wall -Wextra -Wpedantic)
        endif()

        add_executable(orbit_mq_protocol_test tests/orbit_mq_protocol_test.cpp)
        target_link_libraries(orbit_mq_protocol_test PRIVATE orbit_mq_protocol)
        add_test(NAME orbit_mq_protocol_test COMMAND orbit_mq_protocol_test)
        if (MSVC)
            target_compile_options(orbit_mq_protocol_test PRIVATE /W4 /permissive-)
        else()
            target_compile_options(orbit_mq_protocol_test PRIVATE
                -Wall -Wextra -Wpedantic)
        endif()

        add_executable(telemetry_mq_publisher_test tests/telemetry_mq_publisher_test.cpp)
        target_link_libraries(telemetry_mq_publisher_test PRIVATE orbit_prediction_mq_service)
        add_test(NAME telemetry_mq_publisher_test COMMAND telemetry_mq_publisher_test)
        if (MSVC)
            target_compile_options(telemetry_mq_publisher_test PRIVATE /W4 /permissive-)
        else()
            target_compile_options(telemetry_mq_publisher_test PRIVATE
                -Wall -Wextra -Wpedantic)
        endif()

        add_executable(heartbeat_mq_publisher_test tests/heartbeat_mq_publisher_test.cpp)
        target_link_libraries(heartbeat_mq_publisher_test PRIVATE orbit_prediction_mq_service)
        add_test(NAME heartbeat_mq_publisher_test COMMAND heartbeat_mq_publisher_test)
        if (MSVC)
            target_compile_options(heartbeat_mq_publisher_test PRIVATE /W4 /permissive-)
        else()
            target_compile_options(heartbeat_mq_publisher_test PRIVATE
                -Wall -Wextra -Wpedantic)
        endif()
    endif()

    install(TARGETS orbit_prediction_server orbit_mq_cli RUNTIME DESTINATION bin)

endif()
