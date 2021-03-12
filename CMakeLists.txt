cmake_minimum_required(VERSION 3.15)
CMAKE_POLICY(SET CMP0012 NEW)
CMAKE_POLICY(SET CMP0025 NEW)
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.11.0")
    CMAKE_POLICY(SET CMP0072 NEW)
endif()
#============================================================
project(TrackingTelescope)
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
#=========================================================
# Help function to debug CMake
MACRO (DD in)
    MESSAGE(${in}=[${${in}}])
ENDMACRO(DD)
#=========================================================
# Locate sources and headers for this project
FILE(GLOB sources
        ${PROJECT_SOURCE_DIR}/src/*.cc)

FILE(GLOB headers
        ${PROJECT_SOURCE_DIR}/include/*.h)
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
    #MESSAGE("Found Root in ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARY_DIR} ${ROOT_LIBRARIES}")
    #IF(ROOT_INCLUDE_DIR)
    #MESSAGE("Your current Root version is ${ROOT_VERSION}")
    #ELSE()
    IF(NOT ROOT_INCLUDE_DIR)
        MESSAGE("Root not found! Please install Root and make it available in PATH")
    ENDIF()
ENDIF()
#=========================================================
# Add the executable, and link it to the Geant4/ROOT/CLHEP/ITK libraries
ADD_LIBRARY(TrackingTelescopeLib OBJECT ${sources} ${headers})
ADD_EXECUTABLE(TrackingTelescope ${PROJECT_SOURCE_DIR}/src/TrackingTelescope.cxx
        $<TARGET_OBJECTS:TrackingTelescopeLib> )
TARGET_LINK_LIBRARIES(TrackingTelescope ${ROOT_LIBRARIES})

#=========================================================
# We remove the warning option "shadow", because there are tons of
# such warning related to clhep/g4 system of units.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-shadow  -std=c++14")