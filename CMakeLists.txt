cmake_minimum_required(VERSION 3.3)
project(cif2xyz)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package (Eigen3 REQUIRED)
if(FFTW_FOUND)
	message(STATUS "EIGEN3 found (include: ${EIGEN3_INCLUDE_DIR})")
endif(FFTW_FOUND)

set(HEADER_FILES
		unitcell.h)

set(SOURCE_FILES
		main.cpp
        reader.cpp
		symmetryoperation.cpp
        atomsite.cpp
		cellgeometry.cpp
		utilities.cpp
		supercell.cpp)

include_directories(
		${EIGEN3_INCLUDE_DIR})

add_executable(cif2xyz ${SOURCE_FILES} ${HEADER_FILES})