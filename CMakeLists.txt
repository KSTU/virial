cmake_minimum_required(VERSION 2.8)
project(virial)
#Version number
set(virial_VERSION_MAJOR 0)
set(virial_VERSION_MINOR 1)

#Fing cuda
find_package(CUDA QUIET REQUIRED)
message(${CUDA_VERSION})
message(${CMAKE_CXX_COMPILER})
cuda_add_executable(virial main.cu)
