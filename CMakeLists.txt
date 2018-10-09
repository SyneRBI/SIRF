
#========================================================================
# Author: Kris Thielemans
# Copyright 2016 - 2018 University College London
# Copyright 2016 - 2018 Science Technology Facilities Council
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

include(cmake/SetC++Version.cmake)
UseCXX(11)
#set(CMAKE_CXX_STANDARD 11)

PROJECT(SIRF)

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
# require 2.8.3 to get FOLDER properties support
# require 3.3 for descent FindMatlab.cmake
cmake_minimum_required(VERSION 3.3.0)

# set default build-type to Release
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

####### Set Version number etc
set(VERSION_MAJOR   1 )
set(VERSION_MINOR   1 )
set(VERSION_PATCH   1 )
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)

set(SIRF_VERSION
  ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

option(SIRF_INSTALL_DEPENDENCIES "Install dlls etc" WIN32)
####### CMake path
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# If OSX give the advanced option to use absolute paths for shared libraries
if (APPLE)
  option(SHARED_LIBS_ABS_PATH "Force shared libraries to be installed with absolute paths (as opposed to rpaths)" ON)
  mark_as_advanced( SHARED_LIBS_ABS_PATH )  
  if (SHARED_LIBS_ABS_PATH)
    # Set install_name_dir as the absolute path to install_prefix/lib
    GET_FILENAME_COMPONENT(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib REALPATH)
  endif(SHARED_LIBS_ABS_PATH)
endif(APPLE)

####### External packages

#### we need the boost library from boost.org
set(BOOST_ROOT CACHE PATH "root of Boost")
find_package( Boost 1.36.0 REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

#### optional back-ends
# find Python interpreter. Needed for tests, and best to enforce consistency anyway.
find_package(PythonInterp)
find_package(PythonLibs)
set (BUILD_PYTHON ${PYTHONLIBS_FOUND})
if (BUILD_PYTHON)
# PYTHON_DEST_DIR allows the user to select the install destination of the 
# SIRF python modules. PYTHON_DEST_DIR is a cached variable which can be 
# updated on the GUI.
# If PYTHON_DEST_DIR is not set, we will install in ${CMAKE_INSTALL_PREFIX}/python 
  set(PYTHON_DEST_DIR "" CACHE PATH "Directory of the SIRF Python modules")
  if (PYTHON_DEST_DIR)
   set(PYTHON_DEST "${PYTHON_DEST_DIR}")
  else() 
    set(PYTHON_DEST "${CMAKE_INSTALL_PREFIX}/python")
  endif()
  message(STATUS "Python libraries found")
  message(STATUS "SIRF Python modules will be installed in " ${PYTHON_DEST})
  set(PYTHON_STRATEGY "PYTHONPATH" CACHE STRING "\
    PYTHONPATH: prefix PYTHONPATH \n\
    SETUP_PY:   execute ${PYTHON_EXECUTABLE} setup.py install \n\
    CONDA:      do nothing")
endif(BUILD_PYTHON)


message(STATUS "Attempting to find MATLAB (This might launch MATLAB so might take a while)")
FIND_PACKAGE(Matlab COMPONENTS MAIN_PROGRAM)
set (BUILD_MATLAB ${Matlab_FOUND})

if (BUILD_MATLAB)
  matlab_get_mex_suffix("${Matlab_ROOT_DIR}" MATLAB_MEX_EXT)
# MATLAB_DEST_DIR allows the user to select the install destination of the 
# SIRF matlab modules. MATLAB_DEST_DIR is a cached variable which can be 
# updated on the GUI. 
# If MATLAB_DEST_DIR is not set, we will install in ${CMAKE_INSTALL_PREFIX}/matlab 
  set(MATLAB_DEST_DIR "" CACHE PATH "Directory of the SIRF Matlab libraries")
  if (MATLAB_DEST_DIR)
   set(MATLAB_DEST "${MATLAB_DEST_DIR}")
  else() 
    set(MATLAB_DEST "${CMAKE_INSTALL_PREFIX}/matlab")
  endif()
  message(STATUS "Matlab libraries found")
  message(STATUS "SIRF Matlab libraries will be installed in " ${MATLAB_DEST})
endif()

if (NOT BUILD_MATLAB AND NOT BUILD_PYTHON)
  message(FATAL_ERROR "SIRF needs either MATLAB or Python")
endif()

ENABLE_TESTING()

ADD_SUBDIRECTORY(src/iUtilities)
ADD_SUBDIRECTORY(src/xSTIR)
ADD_SUBDIRECTORY(src/xGadgetron)
ADD_SUBDIRECTORY(src/common)
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/setup.py.cmake")

# Optional registration
option(USE_SIRFREG "Build the SIRF registration package" OFF)
if (USE_SIRFREG)
  ADD_SUBDIRECTORY(src/Registration)
endif()
