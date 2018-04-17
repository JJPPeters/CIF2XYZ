cmake_minimum_required(VERSION 3.3)
project(XYZ__)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

include_directories(D:/Programming/Cpp/Eigen/eigen3.3.4)

set(SOURCE_FILES main.cpp CIFReader.cpp CIFReader.h UnitCell.cpp UnitCell.h SymmetryOperation.cpp SymmetryOperation.h Utilities.h AtomSite.cpp AtomSite.h CellGeometry.cpp CellGeometry.h)
add_executable(XYZ__ ${SOURCE_FILES})