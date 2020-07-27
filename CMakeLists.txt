cmake_minimum_required(VERSION 3.0.0)
project(ACESLib VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_library(ACESLib SHARED ACESLib.cpp metal_math.cpp)

