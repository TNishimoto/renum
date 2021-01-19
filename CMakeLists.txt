cmake_minimum_required(VERSION 2.8)

if(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall -pthread -Wextra -static")
else(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall -pthread -Wextra")
endif(CPP_STATIC)

#set(CMAKE_CXX_FLAGS_DEBUG  "-g -O0 -DDEBUG -DDEBUG_PRINT")
set(CMAKE_CXX_FLAGS_DEBUG  "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

#set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

#### Modules
set(MYMOD_SOURCE_DIR ${CMAKE_SOURCE_DIR}/module)
set(MYMOD_BINARY_DIR ${CMAKE_BINARY_DIR}/module)
set(MYMOD_DIR 1)

#add_subdirectory(${MYMOD_SOURCE_DIR}/libdivsufsort ${MYMOD_BINARY_DIR}/libdivsufsort)
#include_directories(${MYMOD_BINARY_DIR}/libdivsufsort/include)

if(COMPILE_LIBDIVSUFSORT)
add_subdirectory(${MYMOD_SOURCE_DIR}/libdivsufsort ${MYMOD_BINARY_DIR}/libdivsufsort)
include_directories(${MYMOD_BINARY_DIR}/libdivsufsort/include)
else()
#libdivsufsort directory is here.
#Please modify it appropriately.
LINK_DIRECTORIES(~/lib)
INCLUDE_DIRECTORIES(~/include)
endif()


#SDSL lib and include directory are here
#Please modify it appropriately.
LINK_DIRECTORIES(~/lib)
include_directories(~/include)


#set(GLOB SOURCES src/minimal_substrings/*.cpp)

# STOOL submodule
FILE(GLOB STOOL_SOURCES module/stool/src/*.cpp module/libdivsufsort/sa.cpp)
INCLUDE_DIRECTORIES(module)

FILE(GLOB STNODE_SOURCES src/basic/*.cpp src/rlbwt/range_distinct/*.cpp)


add_library(fn STATIC ${STOOL_SOURCES} ${STNODE_SOURCES})
target_link_libraries(fn sdsl divsufsort divsufsort64)

#add_library(fn2 STATIC ${STOOL_SOURCES} ${STNODE_SOURCES} src/main/library.cpp)
#target_link_libraries(fn2 sdsl divsufsort divsufsort64)


### Debug
#add_executable(beller.out src/main/beller_main.cpp)
#target_link_libraries(beller.out fn)

add_executable(bbo_maximal_repeat.out src/main/bbo_meximal_repeat_main.cpp)
target_link_libraries(bbo_maximal_repeat.out fn)

add_executable(dfs_maximal_repeat.out src/main/dfs_maximal_repeat_main.cpp)
target_link_libraries(dfs_maximal_repeat.out fn)


add_executable(maximal_repeat.out src/main/maximal_repeat_main.cpp)
target_link_libraries(maximal_repeat.out fn)

add_executable(mus.out src/main/mus_main.cpp)
target_link_libraries(mus.out fn)


add_executable(interval_test.out src/debug/lcp_interval_test_main.cpp)
target_link_libraries(interval_test.out fn)


add_executable(convert_bwt_into_int_vector.out src/main/convert_bwt_into_int_vector_main.cpp)
target_link_libraries(convert_bwt_into_int_vector.out fn)

add_executable(text_to_bwt.out src/main/text_to_bwt_main.cpp)
target_link_libraries(text_to_bwt.out divsufsort divsufsort64)

add_executable(print.out src/main/print.cpp)
target_link_libraries(print.out fn)

#add_executable(bwt_test.out src/debug/bwt_test.cpp)
#target_link_libraries(bwt_test.out fn)
