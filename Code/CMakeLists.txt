cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

project(nanojpeg)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX )

add_executable(nanojpeg app.cpp nanojpeg.hpp mio.hpp cached_bitstream.hpp)


add_compile_options(nanojpeg  -mtune=ivybridge  -mcpu=ivybridge)
if(CMAKE_BUILD_TYPE MATCHES "Release")
    add_compile_options(nanojpeg -msse4.1 -Ofast)
endif()

find_package(Stb REQUIRED)
target_include_directories(nanojpeg PRIVATE ${Stb_INCLUDE_DIR})

find_package(libjpeg-turbo CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE libjpeg-turbo::turbojpeg-static)

