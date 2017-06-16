cmake_minimum_required(VERSION 3.7)
project(Diplom)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp Triangle/triangle.cpp EasyLogging/easylogging++.cc constants.h Gauss/gauss.cpp)
add_executable(Diplom ${SOURCE_FILES})