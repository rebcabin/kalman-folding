cmake_minimum_required(VERSION 3.5)
project(kalmanSample)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -Wall -Werror")

set(CMAKE_SHARED_LINKER_FLAGS "-L /usr/local/lib -lcblas -lLAPACKE -llapack")

# link_directories(/usr/local/lib/)

set(SOURCE_FILES kalmanSample.c)

add_executable(kalmanSample ${SOURCE_FILES})
