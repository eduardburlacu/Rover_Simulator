cmake_minimum_required(VERSION 3.22)
project(lander)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(lander
        Tasks_3_4/try.cpp
        spring.cpp)
