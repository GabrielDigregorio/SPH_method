PROJECT(SPH CXX)
cmake_minimum_required(VERSION 2.7)

FILE(GLOB SRCS *.h *.cpp)
ADD_EXECUTABLE(sph ${SRCS})

# - OpenMP --

find_package(OpenMP REQUIRED)

if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


