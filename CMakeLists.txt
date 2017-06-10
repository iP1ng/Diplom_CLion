cmake_minimum_required(VERSION 3.7)
project(Diplom)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp triangle.cpp easylogging++.cc constants.h gauss.cpp gauss.h )
add_executable(Diplom ${SOURCE_FILES})