
#========================================================================
# Author: Kris Thielemans
# Copyright 2016 - 2019 University College London
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

# Set the CMake policy for SWIG
# https://cmake.org/cmake/help/v3.14/policy/CMP0078.html#policy:CMP0078
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13") 
  # policy introduced in CMake 3.13
  cmake_policy(SET CMP0078 OLD)
endif()
# avoid warning about WIN32 no longer defined in CYGWIN
set(CMAKE_LEGACY_CYGWIN_WIN32 0) 

include(cmake/SetC++Version.cmake)
UseCXX(11)
#set(CMAKE_CXX_STANDARD 11)

PROJECT(SIRF)

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
# require 2.8.3 to get FOLDER properties support
# require 3.3 for descent FindMatlab.cmake
# require 3.9 for compatible FindOPENMP.cmake
cmake_minimum_required(VERSION 3.9.0)

# set default build-type to Release
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING "type of build: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

####### Set Version number etc
set(VERSION_MAJOR   2 )
set(VERSION_MINOR   1 )
set(VERSION_PATCH   0 )
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
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  endif(SHARED_LIBS_ABS_PATH)
endif(APPLE)

# If OSX, we might need to assume pthreads. For this, first try to find boost_thread.
# If that fails, then give a bit of help.
if (APPLE)
  find_package(Boost COMPONENTS thread QUIET)
  if (NOT boost_thread_FOUND)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
  endif()
endif()

####### External packages

#### we need the boost library from boost.org
set(BOOST_ROOT CACHE PATH "root of Boost")
find_package(Boost 1.36.0 COMPONENTS system filesystem thread date_time chrono REQUIRED)
# For Visual Studio we have to disable the auto-linking feature of boost
# where just including a boost file automatically adds it to the linker path.
# Although this sounds great, it sadly breaks because of conflicting file-paths when linking etc etc.
# In any case, we need to add the libraries by hand for other systems.
# See http://stackoverflow.com/questions/32252016/cmake-visual-studio-build-looks-for-wrong-library
add_definitions(-DBOOST_ALL_NO_LIB)

#### optional back-ends
option(DISABLE_PYTHON "Disable building SIRF python support" OFF)
if (DISABLE_PYTHON)
  message(STATUS "Python support disabled")
else(DISABLE_PYTHON)
  # find Python interpreter. Needed for tests, and best to enforce consistency anyway.
  find_package(PythonInterp QUIET)
  find_package(PythonLibs)
  if (PYTHONLIBS_FOUND)
    set(BUILD_PYTHON ON)
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

  else(PYTHONLIBS_FOUND)
    message(WARNING "Pythonlibs not found. Best to set PYTHON_EXECUTABLE if you want Python support. You can set DISABLE_PYTHON to ON to silence this warning.")
  endif(PYTHONLIBS_FOUND)

endif(DISABLE_PYTHON)


option(DISABLE_Matlab "Disable building SIRF matlab support" OFF)
if (DISABLE_Matlab)
  message(STATUS "Matlab support disabled")
else(DISABLE_Matlab)

  FIND_PACKAGE(Matlab QUIET COMPONENTS MAIN_PROGRAM)
  if (NOT Matlab_FOUND)
    message(WARNING "MATLAB not found. Set Matlab_ROOT_DIR if you want it. Set DISABLE_Matlab to ON to silence this warning.")
  else()
    set(BUILD_MATLAB ON)
    message(STATUS "Attempting to find MATLAB Mex extension (This might launch MATLAB so might take a while)")
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

endif(DISABLE_Matlab)

ENABLE_TESTING()

ADD_SUBDIRECTORY(src)

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/setup.py.cmake")