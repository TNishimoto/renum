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

if(NOT TARGET libdivsufsort)
include_directories(${CMAKE_BINARY_DIR}/module/libdivsufsort/include)
add_subdirectory(module/libdivsufsort)
endif()

#set(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")


#SDSL lib and include directory are here
#Please modify it appropriately.
find_library(SDSL_LIBRARY_DIR sdsl REQUIRED)
find_path(SDSL_INCLUDE_DIR sdsl REQUIRED)

LINK_DIRECTORIES(${SDSL_LIBRARY_DIR})
include_directories(${SDSL_INCLUDE_DIR})

INCLUDE_DIRECTORIES(module)


add_executable(maximal_repeat.out main/maximal_repeat_main.cpp)
target_link_libraries(maximal_repeat.out sdsl divsufsort divsufsort64)

add_executable(mus.out main/mus_main.cpp)
target_link_libraries(mus.out sdsl divsufsort divsufsort64)

add_executable(text_to_bwt.out main/text_to_bwt_main.cpp)
target_link_libraries(text_to_bwt.out sdsl divsufsort divsufsort64)

add_executable(print.out main/print.cpp)
target_link_libraries(print.out sdsl divsufsort divsufsort64)

add_executable(exsample.out main/library_sample_main.cpp)
target_link_libraries(exsample.out sdsl divsufsort divsufsort64)





#add_executable(bwt_test.out src/debug/bwt_test.cpp)
#target_link_libraries(bwt_test.out fn)

### For Experiments

add_executable(convert_bwt_into_int_vector.out experiments/convert_bwt_into_int_vector_main.cpp)
target_link_libraries(convert_bwt_into_int_vector.out sdsl divsufsort divsufsort64)

add_executable(bbo_maximal_repeat.out experiments/bbo_meximal_repeat_main.cpp)
target_link_libraries(bbo_maximal_repeat.out sdsl divsufsort divsufsort64)

#add_executable(dfs_maximal_repeat.out experiments/dfs_maximal_repeat_main.cpp)
#target_link_libraries(dfs_maximal_repeat.out fn)



#add_executable(interval_test.out src/debug/lcp_interval_test_main.cpp)
#target_link_libraries(interval_test.out fn)

#add_library(fn2 STATIC ${STOOL_SOURCES} ${STNODE_SOURCES} src/main/library.cpp)
#target_link_libraries(fn2 sdsl divsufsort divsufsort64)


### Debug
#add_executable(beller.out src/main/beller_main.cpp)
#target_link_libraries(beller.out fn)

