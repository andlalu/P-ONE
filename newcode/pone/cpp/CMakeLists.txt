cmake_minimum_required(VERSION 3.27)
project(pone)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(extern/pybind11)
find_package(pybind11 CONFIG)

include_directories(wrap)
include_directories(utils)

set(SOURCE_FILES wrap/wrap.cpp utils/rand.cpp)
pybind11_add_module(pone ${SOURCE_FILES})



