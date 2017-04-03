PROJECT(SPH CXX)
cmake_minimum_required(VERSION 2.7)
include(CheckCXXCompilerFlag)

IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF()

IF(MSVC)
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE)
ENDIF(MSVC)

# --------------------------------------------------------------
# Finds MPI (including MSMPI)
# --------------------------------------------------------------
if(WIN32)  
    FIND_PATH(MPI_INCLUDE_PATH NAMES  mpi.h HINTS "$ENV{MSMPI_INC}")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL 8)
        set(MS_MPI_ARCH_DIR x64)
        find_library(MPI_LIBRARIES msmpi PATHS "$ENV{MSMPI_LIB64}")
    else()
        set(MS_MPI_ARCH_DIR x86)
        find_library(MPI_LIBRARIES msmpi PATHS "$ENV{MSMPI_LIB32}")
    endif()
    message(STATUS "MPI_LIBRARIES = ${MPI_LIBRARIES}")
    IF( CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        add_definitions(-DMSMPI_NO_SAL -D"__int64=long long") # MinGW
    ENDIF()
else()
    find_package(MPI REQUIRED)
endif()
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
# --------------------------------------------------------------

check_cxx_compiler_flag("-std=c++11" STD11CHECK)
IF(STD11CHECK)
    ADD_DEFINITIONS(-std=c++11 )
ENDIF()

FILE(GLOB SRCS1  CPP_Interface/*.cpp)
#FILE(GLOB SRCS1 CPP_Interface/ParaView.cpp CPP_Interface/writeField.cpp)
FILE(GLOB SRCS1 CPP_Interface/*.h CPP_Interface/*.cpp)
FILE(GLOB SRCS2 CPP_Physics/*.h CPP_Physics/*.cpp)
FILE(GLOB SRCS3 CPP_Tools/*.h CPP_Tools/*.cpp)
FILE(GLOB SRCS4 Headers/*.h Headers/*.cpp)

INCLUDE_DIRECTORIES(CPP_Interface CPP_Physics CPP_Tools Headers)

ADD_EXECUTABLE(sph ${SRCS1} ${SRCS2} ${SRCS3} ${SRCS4} CPP_Main/SPH.cpp)
#ADD_EXECUTABLE(neighbors ${SRCS1} ${SRCS2} ${SRCS3} ${SRCS4} CPP_Main/Neighborhood_performance.cpp)

target_link_libraries(sph ${MPI_LIBRARIES})

IF(MINGW)
    TARGET_LINK_LIBRARIES(sph psapi) # for "GetProcessMemoryInfo"
    TARGET_LINK_LIBRARIES(neighbors psapi)
ENDIF(MINGW)

# - OpenMP --
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 OLD)  # "if() recognizes numbers and boolean constants."
endif()
find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
