cmake_minimum_required(VERSION 3.3)
project(cif2xyz)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

find_package (Eigen3 REQUIRED)
if(FFTW_FOUND)
	message(STATUS "EIGEN3 found (include: ${EIGEN3_INCLUDE_DIR})")
endif(FFTW_FOUND)

set(HEADER_FILES
		Utilities.h)

set(SOURCE_FILES
		main.cpp
		CIFReader.cpp
		UnitCell.cpp
		SymmetryOperation.cpp
		AtomSite.cpp
		CellGeometry.cpp)

include_directories(
		${EIGEN3_INCLUDE_DIR})

add_executable(cif2xyz ${SOURCE_FILES} ${HEADER_FILES})