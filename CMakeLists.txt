cmake_minimum_required(VERSION 3.30)
set(CMAKE_CXX_STANDARD 17)
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang|GNU")
    set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
project(nanojpeg)

#set(CMAKE_NASM_LINK_EXECUTABLE "ld <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES>")
string(TOUPPER ${CMAKE_BUILD_TYPE} CONFIG)

set(CMAKE_ASM_NASM_OBJECT_FORMAT win64)
set(CMAKE_ASM_NASM_FLAGS_DEBUG_INIT "-gcv8")
set(CMAKE_ASM_NASM_FLAGS_RELWITHDEBINFO_INIT "-gcv8")
enable_language(ASM_NASM)
set(CMAKE_ASM_NASM_FLAGS "${CMAKE_ASM_NASM_FLAGS} -DWIN64")


if(CMAKE_CXX_SIMULATE_ID STREQUAL "")
    set(COMPILER_INFO "${CMAKE_CXX_COMPILER_ID}")
else()
    set(COMPILER_INFO "${CMAKE_CXX_SIMULATE_ID}")
endif()
add_executable(nanojpeg Code/app.cpp Code/idct.asm Code/nanojpeg.hpp Code/profiling.hpp Code/converter.hpp)
set_source_files_properties(Code/idct.asm PROPERTIES LANGUAGE ASM_NASM)
if(MSVC)
target_compile_definitions(nanojpeg PRIVATE $<$<COMPILE_LANGUAGE:CXX>:_USE_MATH_DEFINES NOMINMAX>)
endif()



if(COMPILER_INFO MATCHES "Clang|GNU")
    #string(APPEND  CMAKE_CXX_FLAGS " -mssse3 -msse4.1 -mavx -march=sandybridge -O3 -ffast-math ")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mssse3 -msse4.1 -mavx -march=sandybridge -mtune=sandybridge>)
    if( CONFIG STREQUAL "RELEASE" )
        add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-O3 -ffast-math>)
    endif()
    if(COMPILER_INFO MATCHES "Clang")
        #add_compile_options( -flto=thin)
        #string(APPEND  CMAKE_CXX_FLAGS  " -flto=thin ")
        #add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-flto=thin>)
    endif()
    if(COMPILER_INFO MATCHES "GNU")
        #add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-flto>)
    endif()

endif()

if(COMPILER_INFO MATCHES "MSVC")
    #string(APPEND  CMAKE_CXX_FLAGS  " -mssse3 -msse4.1 -mavx ")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mssse3 -msse4.1 -mavx>)
endif()


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

find_package(benchmark CONFIG)
if(benchmark)
add_executable(bench benchmark/app.cpp)
target_link_libraries(bench PRIVATE benchmark::benchmark benchmark::benchmark_main)
endif()