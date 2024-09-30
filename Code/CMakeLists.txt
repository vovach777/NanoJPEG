cmake_minimum_required(VERSION 3.15)
project(nanojpeg)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX )

add_executable(nanojpeg app.cpp nanojpeg.hpp)
set_target_properties(nanojpeg PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    target_compile_options(nanojpeg PRIVATE -O3 -flto  -march=native -mavx  -fopenmp-simd)
endif()

find_package(Stb REQUIRED)
target_include_directories(nanojpeg PRIVATE ${Stb_INCLUDE_DIR})

find_package(libjpeg-turbo CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE libjpeg-turbo::turbojpeg-static)

find_package(mio CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE mio::mio)

option(USE_XSIMD "Build with xsimd support" OFF)

if (USE_XSIMD)
    find_package(xsimd CONFIG REQUIRED)
    if (xsimd_FOUND)
        target_link_libraries(nanojpeg PRIVATE xsimd)
        target_compile_definitions(nanojpeg PRIVATE HAS_XSIMD)
    endif()
endif()

find_path(SIMDE_INCLUDE_DIRS "simde/x86/avx.h")
target_include_directories(nanojpeg PRIVATE ${SIMDE_INCLUDE_DIRS})