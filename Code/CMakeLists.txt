cmake_minimum_required(VERSION 3.30)
project(nanojpeg)
#set(CMAKE_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})
set(CMAKE_CXX_STANDARD 17)
string(TOUPPER ${CMAKE_BUILD_TYPE} CONFIG)
add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX )

add_executable(nanojpeg app.cpp nanojpeg.hpp)
set_target_properties(nanojpeg PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

string(APPEND CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -msse4.2 -mavx  -fopenmp-simd")

string(APPEND CMAKE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")

find_package(Stb REQUIRED)
target_include_directories(nanojpeg PRIVATE ${Stb_INCLUDE_DIR})

find_package(libjpeg-turbo CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE libjpeg-turbo::turbojpeg-static)

find_package(mio CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE mio::mio)
message(STATUS "c++ opt_${CONFIG}:  ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CONFIG}}")
message(STATUS "link opt_${CONFIG}:  ${CMAKE_EXE_LINKER_FLAGS_INIT}")