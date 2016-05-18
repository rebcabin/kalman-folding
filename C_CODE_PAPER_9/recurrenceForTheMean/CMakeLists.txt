cmake_minimum_required(VERSION 3.5)
project(recurrenceForTheMean)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Werror")

set(SOURCE_FILES main.c)
add_executable(recurrenceForTheMean ${SOURCE_FILES})