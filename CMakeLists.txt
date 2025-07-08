cmake_minimum_required(VERSION 2.8)

#Project Information
project(renum)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


if(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-m64 -Wall -pthread -Wextra -static")
else(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-m64 -Wall -pthread -Wextra")
endif(CPP_STATIC)

#set(CMAKE_CXX_FLAGS_DEBUG  "-g -O0 -DDEBUG -DDEBUG_PRINT")
set(CMAKE_CXX_FLAGS_DEBUG  "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

#set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

#### Modules
set(MYMOD_SOURCE_DIR ${CMAKE_SOURCE_DIR}/module)
set(MYMOD_BINARY_DIR ${CMAKE_BINARY_DIR}/module)
set(MYMOD_DIR 1)

#SDSL lib and include directory are here
#Please modify it appropriately.
find_library(SDSL_LIBRARY_DIR sdsl REQUIRED)
find_path(SDSL_INCLUDE_DIR sdsl REQUIRED)

LINK_DIRECTORIES(${SDSL_LIBRARY_DIR})
include_directories(${SDSL_INCLUDE_DIR})
#include_directories(~/include)

if(NOT TARGET stool)
add_subdirectory(module/stool)
endif()



#set(GLOB SOURCES src/minimal_substrings/*.cpp)

# STOOL submodule
FILE(GLOB STOOL_SOURCES module/libdivsufsort/sa.cpp)
INCLUDE_DIRECTORIES(module)
include_directories(${PROJECT_SOURCE_DIR}/include)


FILE(GLOB STNODE_SOURCES src/lib/basic/*.cpp src/lib/rlbwt/*.cpp)


add_library(fn STATIC ${STOOL_SOURCES} ${STNODE_SOURCES})
target_link_libraries(fn stool sdsl divsufsort divsufsort64)



add_executable(maximal_repeat.out src/main/maximal_repeat_main.cpp)
target_link_libraries(maximal_repeat.out fn)

add_executable(mus.out src/main/mus_main.cpp)
target_link_libraries(mus.out fn)

add_executable(text_to_bwt.out src/main/text_to_bwt_main.cpp)
target_link_libraries(text_to_bwt.out fn)

add_executable(print.out src/main/print.cpp)
target_link_libraries(print.out fn)

add_executable(exsample.out src/main/library_sample_main.cpp)
target_link_libraries(exsample.out fn)





#add_executable(bwt_test.out src/debug/bwt_test.cpp)
#target_link_libraries(bwt_test.out fn)

### For Experiments

add_executable(convert_bwt_into_int_vector.out src/experiments/convert_bwt_into_int_vector_main.cpp)
target_link_libraries(convert_bwt_into_int_vector.out fn)

add_executable(bbo_maximal_repeat.out src/experiments/bbo_meximal_repeat_main.cpp)
target_link_libraries(bbo_maximal_repeat.out fn)

#add_executable(dfs_maximal_repeat.out src/experiments/dfs_maximal_repeat_main.cpp)
#target_link_libraries(dfs_maximal_repeat.out fn)



#add_executable(interval_test.out src/debug/lcp_interval_test_main.cpp)
#target_link_libraries(interval_test.out fn)

#add_library(fn2 STATIC ${STOOL_SOURCES} ${STNODE_SOURCES} src/main/library.cpp)
#target_link_libraries(fn2 sdsl divsufsort divsufsort64)


### Debug
#add_executable(beller.out src/main/beller_main.cpp)
#target_link_libraries(beller.out fn)

