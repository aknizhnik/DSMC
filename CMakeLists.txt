cmake_minimum_required(VERSION 3.19)
project(DSMC)

set(CMAKE_CXX_STANDARD 17)
include_directories(/opt/local/include/eigen3.2)

add_executable(DSMC main.cpp dsmc.cpp)