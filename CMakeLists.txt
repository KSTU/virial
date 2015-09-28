cmake_minimum_required(VERSION 2.8)
project(virial)
#Version number
set(virial_VERSION_MAJOR 0)
set(virial_VERSION_MINOR 1)

#Fing cuda
find_package(CUDA QUIET REQUIRED)
message("cuda version" ${CUDA_VERSION})
message("cmake version" ${CMAKE_CXX_COMPILER})
cuda_add_executable(virial main.cu)
cuda_add_library(messages STATIC messages.cu)
cuda_add_library(input STATIC input.cu)
cuda_add_library(integ STATIC integrate.cu)
target_link_libraries(virial messages input integ)
