PROJECT(SPH CXX)
cmake_minimum_required(VERSION 2.7)
IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF()

IF(MSVC)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE)
ENDIF(MSVC)
IF((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
    ADD_DEFINITIONS(-std=c++11 )
ENDIF()

#FILE(GLOB SRCS1 CPP_Interface/*.h CPP_Interface/*.cpp)
#FILE(GLOB SRCS1 CPP_Interface/ParaView.cpp CPP_Interface/writeField.cpp)
FILE(GLOB SRCS1 CPP_Interface/writeField.cpp)
FILE(GLOB SRCS2 CPP_Physics/*.h CPP_Physics/*.cpp)
FILE(GLOB SRCS3 CPP_Tools/*.h CPP_Tools/*.cpp)
FILE(GLOB SRCS4 Headers/*.h Headers/*.cpp)

INCLUDE_DIRECTORIES(CPP_Interface CPP_Physics CPP_Tools Headers)

ADD_EXECUTABLE(sph ${SRCS1} ${SRCS2} ${SRCS3} ${SRCS4} CPP_Main/SPH.cpp)
ADD_EXECUTABLE(neighbors ${SRCS1} ${SRCS2} ${SRCS3} ${SRCS4} CPP_Main/Neighborhood_performance.cpp)



# - OpenMP --
find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
