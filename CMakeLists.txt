cmake_minimum_required(VERSION 3.15)
project(project)

set(CMAKE_CXX_STANDARD 17)

add_executable(project main.cpp matrix.cpp matrix.h)
#add_executable(project main.cpp matrix.cpp matrix.h k_means.h k_means.cpp)
#add_executable(project main.cpp csv_reader.cpp csv_reader.h)