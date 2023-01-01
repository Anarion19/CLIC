cmake_minimum_required(VERSION 3.14)
project(CLIC)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O4 -march=native -fext-numeric-literals -std=c++17 -lboost_system -fopenmp")

include_directories(.)

add_executable(CLIC
        data_holders.h
        main.cpp
        new_case.h
        simulation.h
        transfer.h
        transformation.h genetic.h)

target_link_libraries (CLIC Eigen3::Eigen)