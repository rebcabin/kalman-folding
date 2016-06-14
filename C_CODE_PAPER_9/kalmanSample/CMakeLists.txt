cmake_minimum_required(VERSION 3.5)

project(kalmanSample)

set(SOURCE_FILES kalmanSample.c)

message (STATUS "Source files: " ${SOURCE_FILES})

add_executable(kalmanSample ${SOURCE_FILES})

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -Wall -Werror")

find_library (BLAS_LIB    NAMES blas     PATHS "/usr/local/lib")
find_library (CBLAS_LIB   NAMES cblas    PATHS "/usr/local/lib")
find_library (LAPACKE_LIB NAMES lapacke  PATHS "/usr/local/lib")
find_library (LAPACK_LIB  NAMES lapack   PATHS "/usr/local/lib")
find_library (GFORT_LIB   NAMES gfortran PATHS "/usr/local/Cellar/gcc/5.3.0/lib/gcc/5")

target_link_libraries (kalmanSample
  ${BLAS_LIB}
  ${CBLAS_LIB}
  ${LAPACKE_LIB}
  ${LAPACK_LIB}
  ${GFORT_LIB}
  )

