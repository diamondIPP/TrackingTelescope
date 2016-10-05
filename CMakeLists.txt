cmake_minimum_required(VERSION 3.3)
project(TrackingTelescope)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES bin/TrackingTelescope.cc)
add_executable(TrackingTelescope ${SOURCE_FILES})

include_directories(interface)
add_subdirectory(src)
#ADD_CUSTOM_TARGET(bla COMMAND make -C Makefile_ubuntu)