cmake_minimum_required(VERSION 3.27)
project(pone)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED [COMPONENTS random])
target_link_libraries(main PRIVATE Boost::boost Boost::random)

add_library(pone SHARED pone.cpp)
