cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(simPol)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O2 -fopenmp")

find_package(OpenMP REQUIRED)

add_subdirectory(HighFive)
add_executable(simPol simPol.cpp)
target_link_libraries(simPol PUBLIC OpenMP::OpenMP_CXX HighFive)