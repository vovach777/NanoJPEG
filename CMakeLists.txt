cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 17)
project(nanojpeg)
add_compile_definitions(_USE_MATH_DEFINES NOMINMAX )
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang|IntelLLVM")
    set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

string(TOUPPER ${CMAKE_BUILD_TYPE} CONFIG)

if(CMAKE_CXX_SIMULATE_ID STREQUAL "")
    set(COMPILER_INFO "${CMAKE_CXX_COMPILER_ID}")
else()
    set(COMPILER_INFO "${CMAKE_CXX_SIMULATE_ID}")
endif()
if(COMPILER_INFO MATCHES "Clang|GNU|IntelLLVM")
    add_compile_options(-mssse3 -msse4.1 -march=core2 -mtune=core2 $<$<CONFIG:Release>-O3 -ffast-math>)
    if(COMPILER_INFO MATCHES "Clang")
        add_compile_options( -flto=thin)
    else()
        #add_compile_options( -flto )
    endif()
endif()
if(COMPILER_INFO MATCHES "MSVC")
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|IntelLLVM")
        add_compile_options( -mssse3 -msse4.1)
        if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
        add_compile_options( $<$<CONFIG:Release>:/fast> )
        endif()
    endif()

endif()
add_executable(nanojpeg Code/app.cpp Code/nanojpeg.hpp Code/profiling.hpp Code/converter.hpp)

find_package(Stb REQUIRED)
target_include_directories(nanojpeg PRIVATE ${Stb_INCLUDE_DIR})

find_package(libjpeg-turbo CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE libjpeg-turbo::turbojpeg-static)

find_package(mio CONFIG REQUIRED)
target_link_libraries(nanojpeg PRIVATE mio::mio)
# string(REPLACE "-mavx" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# string(REPLACE "-msse4.1" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# string(REPLACE "-march=native" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

message(STATUS "${CMAKE_CXX_COMPILER_ID} (${CMAKE_CXX_SIMULATE_ID})")
message(STATUS "c++ opt_${CONFIG}:  [${CMAKE_CXX_FLAGS}] [${CMAKE_CXX_FLAGS_${CONFIG}}]")
message(STATUS "link opt_${CONFIG}:  [${CMAKE_EXE_LINKER_FLAGS_RELEASE}]")


string(APPEND  CMAKE_EXE_LINKER_FLAGS_RELEASE "-v")

find_package(benchmark CONFIG REQUIRED)
add_executable(bench benchmark/app.cpp)
target_link_libraries(bench PRIVATE benchmark::benchmark benchmark::benchmark_main)

target_compile_options(bench PRIVATE -mavx)
