CMAKE_MINIMUM_REQUIRED(VERSION 3.15)
#============================================================
project(TrackingTelescope)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes")
#=========================================================
# Locate sources and headers for this project
FILE(GLOB sources  ${PROJECT_SOURCE_DIR}/src/*.cc)
FILE(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
#=========================================================
# Setup include directory for this project
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
#=========================================================
# Option for Root use
SET(USE_ROOT ON)
#=========================================================
# ROOT
IF(USE_ROOT)
    FIND_PACKAGE(ROOT REQUIRED COMPONENTS TSpectrum)
    INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
    SET(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
    IF(NOT ROOT_INCLUDE_DIR)
        MESSAGE(FATAL_ERROR "Root not found! Please install Root and make it available in PATH")
    ENDIF()
ENDIF()
#=========================================================
# Add the executable, and link it to the ROOT libraries
ADD_LIBRARY(TrackingTelescopeLib OBJECT ${sources} ${headers})
ADD_EXECUTABLE(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/src/TrackingTelescope.cxx $<TARGET_OBJECTS:TrackingTelescopeLib> )
SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")  # put exe to project dir
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${ROOT_LIBRARIES})
