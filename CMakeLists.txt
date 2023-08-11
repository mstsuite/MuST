
cmake_minimum_required(VERSION 3.18)

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
message(STATUS "CMake ${CMAKE_VERSION} with max. ${Ncpu} threads")

project(MuST LANGUAGES Fortran C)

option(BUILD_MST "Build MST Code" ON)
message(STATUS "Build MST Code: " ${BUILD_MST})
add_subdirectory(MST)