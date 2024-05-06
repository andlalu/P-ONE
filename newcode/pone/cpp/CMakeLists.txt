cmake_minimum_required(VERSION 3.27)
project(pone)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(extern/boostlocal)
find_package(Boost COMPONENTS random)

add_subdirectory("utils")
include_directories("utils")

add_subdirectory(extern/pybind11)
find_package(pybind11 CONFIG)

add_subdirectory("wrap")
include_directories("wrap")
pybind11_add_module(pone wrap/wrap.cpp)



