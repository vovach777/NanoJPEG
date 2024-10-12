cmake_minimum_required(VERSION 3.20)
set(CMAKE_CXX_STANDARD 17)
string(TOUPPER ${CMAKE_BUILD_TYPE} CONFIG)
add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX )
project(nanojpeg)
add_compile_options( -mavx2 -march=native)
add_executable(nanojpeg Code/app.cpp Code/nanojpeg.hpp Code/profiling.hpp Code/converter.hpp)

find_package(Stb REQUIRED)
target_include_directories(nanojpeg PRIVATE ${Stb_INCLUDE_DIR})

find_package(libjpeg-turbo CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE libjpeg-turbo::turbojpeg-static)

find_package(mio CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE mio::mio)

message(STATUS "c++ opt_${CONFIG}:  [${CMAKE_CXX_FLAGS}] [${CMAKE_CXX_FLAGS_${CONFIG}}]")
message(STATUS "link opt_${CONFIG}:  ${CMAKE_EXE_LINKER_FLAGS_INIT}")