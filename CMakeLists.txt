cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(simPol)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -g -fopenmp -DHIGHFIVE_USE_BOOST=OFF")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O2 -fopenmp -DHIGHFIVE_USE_BOOST=OFF")

find_package(OpenMP REQUIRED)

find_package(HighFive REQUIRED)
add_executable(simPol simPol.cpp)
target_link_libraries(simPol PUBLIC OpenMP::OpenMP_CXX HighFive)