# This is a CMake example for Python
cmake_minimum_required(VERSION 3.8)

project(simpleExample)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,/usr/local/lib -lm -ladept")


if(NOT CMAKE_BUILD_TYPE)
 set(CMAKE_BUILD_TYPE Release)
endif()
#Set the cmake flags as O3 always
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -O3 -DADEPT_RECORDING_PAUSABLE")

FIND_PACKAGE(SWIG REQUIRED)
FIND_PACKAGE(PythonLibs)


INCLUDE(${SWIG_USE_FILE})
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES("/usr/local/include")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}")


LINK_DIRECTORIES("/usr/local/lib")


SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(simpleExample.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_LIBRARY(simpleExample LANGUAGE python SOURCES  solveODE.cpp simpleExample.i)
SWIG_LINK_LIBRARIES(simpleExample adept m ${PYTHON_LIBRARIES})