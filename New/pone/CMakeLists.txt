cmake_minimum_required(VERSION 3.27)
project(pone)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost COMPONENTS random)

add_library(pone SHARED pone.cpp)

target_link_libraries(pone Boost::random)


