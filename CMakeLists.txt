cmake_minimum_required(VERSION 2.8)


#if(APPLE)
#set(CMAKE_CXX_FLAGS "-std=c++11  -m64 -w")
#else(APPLE)
#set(CMAKE_CXX_FLAGS "-std=c++11  -m64 -w -static")
#endif(APPLE)

if(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall -pthread -Wextra -static")
else(CPP_STATIC)
set(CMAKE_CXX_FLAGS "-std=c++11 -m64 -Wall -pthread -Wextra")
endif(CPP_STATIC)

#set(CMAKE_CXX_FLAGS_DEBUG  "-g -O0 -DDEBUG -DDEBUG_PRINT")
set(CMAKE_CXX_FLAGS_DEBUG  "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

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

#add_subdirectory(${MYMOD_SOURCE_DIR}/Basics ${MYMOD_BINARY_DIR}/Basics)
#add_subdirectory(${MYMOD_SOURCE_DIR}/BTree ${MYMOD_BINARY_DIR}/BTree)
#add_subdirectory(${MYMOD_SOURCE_DIR}/OnlineRlbwt ${MYMOD_BINARY_DIR}/OnlineRlbwt)

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

#add_library(fn2 STATIC ${STOOL_SOURCES})
#target_link_libraries(fn sdsl divsufsort divsufsort64)


### Main executive files
#add_executable(enumSubstring.out src/main/enumSubstring.cpp)
#target_link_libraries(enumSubstring.out fn)
#
#add_executable(enumMaximalSubstring.out src/main/enumMaximalSubstring.cpp)
#target_link_libraries(enumMaximalSubstring.out fn)
#add_executable(enumLineMaximalSubstring.out src/main/enumLineMaximalSubstring.cpp )
#target_link_libraries(enumLineMaximalSubstring.out fn)
#add_executable(enumMaximalSubstringWithRLBWT.out src/main/enumMaximalSubstringWithRLBWT.cpp )
#target_link_libraries(enumMaximalSubstringWithRLBWT.out fn)
#
#add_executable(print.out src/main/print.cpp)
#target_link_libraries(print.out fn)
#
#
#add_executable(lisInLcpArray.out src/debug/lis_in_lcp_array.cpp )
#target_link_libraries(lisInLcpArray.out fn)
#add_executable(online_rlbwt.out ${MYMOD_SOURCE_DIR}/rlbwt_iterator/src/main/compress_main.cpp)
#target_link_libraries(online_rlbwt.out fn)


### Debug
#add_executable(beller.out src/main/beller_main.cpp)
#target_link_libraries(beller.out fn)

add_executable(hyper.out src/main/hyper_main.cpp)
target_link_libraries(hyper.out fn)



#add_executable(rlesa.out src/debug/rlesa_main.cpp )
#target_link_libraries(rlesa.out fn)

#add_executable(enumMinimalSubstring.out src/main/enumMinimalSubstring.cpp )
#target_link_libraries(enumMinimalSubstring.out fn)

#add_executable(enumKaiMaximalSubstring.out src/main/enumKaiMaximalSubstring.cpp)
#target_link_libraries(enumKaiMaximalSubstring.out fn)




#add_executable(test.out src/debug/test_main.cpp )
#target_link_libraries(test.out fn)

#add_executable(sdsl.out src/debug/sdsl_main.cpp)
#target_link_libraries(sdsl.out fn)

#add_executable(interval_test.out src/debug/interval_test_main.cpp )
#target_link_libraries(interval_test.out fn)

#add_executable(depth.out src/debug/depth_main.cpp )
#target_link_libraries(depth.out fn)

#add_executable(array_test.out src/debug/array_test_main.cpp )
#target_link_libraries(array_test.out fn)

#add_executable(lcp_interval_test.out src/debug/lcp_interval_test_main.cpp )
#target_link_libraries(lcp_interval_test.out fn)



