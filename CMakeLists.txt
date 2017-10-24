cmake_minimum_required(VERSION 3.9)
project(Gauss)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp Matrix.h FieldNumber.h RealNumber.h)
add_executable(Gauss ${SOURCE_FILES})