#========================================================================
# Author: Kris Thielemans
# Copyright 2016 University College London
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#=========================================================================

# cmake file for building SIRF. See the SIRF User's Guide and http://www.cmake.org.

# avoid warning about WIN32 no longer defined in CYGWIN
set(CMAKE_LEGACY_CYGWIN_WIN32 0) 

PROJECT(SIRF)

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
# require 2.8.3 to get FOLDER properties support
cmake_minimum_required(VERSION 2.8.3)

# add project source to cmake path such that it can use our find_package modules and .cmake files
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/src/cmake)

# set default build-type to Release
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

####### Set Version number etc
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   1   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   0   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

set(SIRF_VERSION
  ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

####### CMake path
set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

####### External packages

#### we need the boost library from boost.org
set(BOOST_ROOT CACHE PATH "root of Boost")
find_package( Boost 1.36.0 REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

#### optional back-ends
FIND_PACKAGE(PythonLibs REQUIRED)
set (BUILD_PYTHON ${PythonLibs_FOUND})
set(PYTHON_DEST ${CMAKE_INSTALL_PREFIX}/python CACHE PATH "Destination for python module")

FIND_PACKAGE(MATLAB)
set (BUILD_MATLAB ${MATLAB_FOUND})
set(MATLAB_DEST ${CMAKE_INSTALL_PREFIX}/matlab CACHE PATH "Destination for matlab module")

#### enable support for ctest
ENABLE_TESTING()

ADD_SUBDIRECTORY(iUtilities)
ADD_SUBDIRECTORY(xSTIR)
ADD_SUBDIRECTORY(xGadgetron)
