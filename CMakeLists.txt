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


add_library(fn STATIC ${STOOL_SOURCES})
target_link_libraries(fn sdsl divsufsort divsufsort64)

### Debug
add_executable(beller.out src/main/beller_main.cpp)
target_link_libraries(beller.out fn)

add_executable(beller2.out src/main/beller2_main.cpp)
target_link_libraries(beller2.out fn)


add_executable(hyper.out src/main/hyper_main.cpp)
target_link_libraries(hyper.out fn)

#add_executable(elias_test.out module/stool/src/test/elias_fano_test_main.cpp)
#target_link_libraries(elias_test.out fn)

#add_executable(ptest.out src/main/ptest_main.cpp)
#target_link_libraries(ptest.out fn)

#add_executable(deque.out src/debug/bit_deque_main.cpp)
#target_link_libraries(deque.out fn)

add_executable(interval_test.out src/debug/lcp_interval_test_main.cpp)
target_link_libraries(interval_test.out fn)

add_executable(maximal_repeat_test.out src/debug/maximal_repeat_test_main.cpp)
target_link_libraries(maximal_repeat_test.out fn)
