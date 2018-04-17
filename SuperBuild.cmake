#========================================================================
# Author: Benjamin A Thomas
# Author: Edoardo Pasca
# Copyright 2017 University College London
# Copyright 2017 Science Technology Facilities Council
#
# This file is part of the CCP PETMR Synergistic Image Reconstruction Framework (SIRF) SuperBuild.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#=========================================================================
if (WIN32)
 # used to check for CMAKE_GENERATOR_PLATFORM but that no longer works in 3.10
 if(NOT "x_${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x_x64")
    message(STATUS "CMAKE_GENERATOR: ${CMAKE_GENERATOR}")
    message(STATUS "CMAKE_GENERATOR_PLATFORM: ${CMAKE_GENERATOR_PLATFORM}")
    message(STATUS "CMAKE_VS_PLATFORM_NAME: ${CMAKE_VS_PLATFORM_NAME}")
    message( FATAL_ERROR "The SuperBuild currently has Win64 hard-wired for dependent libraries. Please use a Win64 generator/toolset. Currently using platform '${CMAKE_VS_PLATFORM_NAME}'.")
 endif()
endif()

set(SUPERBUILD_WORK_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH
    "The path for downloading external source directories" )
set(SOURCE_ROOT_DIR "${CMAKE_CURRENT_BINARY_DIR}/sources/"  CACHE PATH "blabla")

mark_as_advanced( SOURCE_DOWNLOAD_CACHE )

set(externalProjName ${PRIMARY_PROJECT_NAME})
set(proj ${PRIMARY_PROJECT_NAME})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/INSTALL" CACHE PATH "Prefix for path for installation" FORCE)
endif()

# If OSX give the advanced option to use absolute paths for shared libraries
if (APPLE)
  option(SHARED_LIBS_ABS_PATH "Force shared libraries to be installed with absolute paths (as opposed to rpaths)" ON)
  mark_as_advanced( SHARED_LIBS_ABS_PATH )
  if (SHARED_LIBS_ABS_PATH)
    # Set install_name_dir as the absolute path to install_prefix/lib
    GET_FILENAME_COMPONENT(CMAKE_INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/lib REALPATH)
    # Mark it as superbuild so that it gets passed to all dependencies
    mark_as_superbuild( PROJECTS ALL_PROJECTS VARS CMAKE_INSTALL_NAME_DIR:PATH )
  endif(SHARED_LIBS_ABS_PATH)
endif(APPLE)

set (SUPERBUILD_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

include(ExternalProject)

set(EXTERNAL_PROJECT_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE INTERNAL "Default build type for support libraries")
message(STATUS "EXTERNAL_PROJECT_BUILD_TYPE: ${EXTERNAL_PROJECT_BUILD_TYPE}")

# Make sure that some CMake variables are passed to all dependencies
mark_as_superbuild(
   PROJECTS ALL_PROJECTS
   VARS CMAKE_GENERATOR:STRING CMAKE_GENERATOR_PLATFORM:STRING CMAKE_GENERATOR_TOOLSET:STRING
        CMAKE_C_COMPILER:FILEPATH CMAKE_CXX_COMPILER:FILEPATH
        CMAKE_INSTALL_PREFIX:PATH
)

# Attempt to make Python settings consistent
set(PYVER 0 CACHE STRING "Python version")
if(PYVER EQUAL 0)
  find_package(PythonInterp)
else()
  find_package(PythonInterp ${PYVER})
endif()
if (PYTHONINTERP_FOUND)
  set(Python_ADDITIONAL_VERSIONS ${PYTHON_VERSION_STRING})
  message(STATUS "Found PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")
  message(STATUS "Python version ${PYTHON_VERSION_STRING}")
endif()
# find_package(PythonLibs ${PYTHON_VERSION_STRING})
if(PYVER EQUAL 0)
  find_package(PythonLibs)
else()
  find_package(PythonLibs ${PYVER})
endif()
if (PYTHONLIBS_FOUND)
  message(STATUS "Found PYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIRS}")
  message(STATUS "Found PYTHON_LIBRARIES=${PYTHON_LIBRARIES}")
endif()

# Find Matlab
set(Matlab_ROOT_DIR $ENV{Matlab_ROOT_DIR} CACHE PATH "Path to Matlab root directory" )
# Note that we need the main program for the configuration files and the tests)
find_package(Matlab COMPONENTS MAIN_PROGRAM)

if (UNIX AND NOT APPLE)
  option(USE_SYSTEM_Boost "Build using an external version of Boost" OFF)
else()
  option(USE_SYSTEM_Boost "Build using an external version of Boost" ON)
endif()
option(USE_SYSTEM_STIR "Build using an external version of STIR" OFF)
option(USE_SYSTEM_HDF5 "Build using an external version of HDF5" OFF)
option(USE_SYSTEM_ISMRMRD "Build using an external version of ISMRMRD" OFF)
option(USE_SYSTEM_siemens_to_ismrmrd "Build using an external version of siemens_to_ismrmrd" OFF)
option(USE_SYSTEM_FFTW3 "Build using an external version of fftw" OFF)
option(USE_SYSTEM_Armadillo "Build using an external version of Armadillo" OFF)
option(USE_SYSTEM_SWIG "Build using an external version of SWIG" OFF)
#option(USE_SYSTEM_Gadgetron "Build using an external version of Gadgetron" OFF)
option(USE_SYSTEM_SIRF "Build using an external version of SIRF" OFF)
option(USE_SYSTEM_GTest "Build using an external version of GTest" OFF)
option(BUILD_STIR_WITH_OPENMP "Build STIR with OpenMP acceleration" OFF)

if (WIN32)
  set(build_Gadgetron_default OFF)
else()
  set(build_Gadgetron_default ON)
endif()

option(BUILD_GADGETRON "Build Gadgetron" ${build_Gadgetron_default})
option(BUILD_siemens_to_ismrmrd "Build siemens_to_ismrmrd" OFF)
option(BUILD_petmr_rd_tools "Build petmr_rd_tools" OFF)

if (BUILD_petmr_rd_tools)
    set(USE_ITK ON CACHE BOOL "Use ITK" FORCE)
    option(USE_SYSTEM_glog "Build using an external version of glog" OFF)
endif()

# ITK
option(USE_ITK "Use ITK" OFF)
if (USE_ITK)
  option(USE_SYSTEM_ITK "Build using an external version of ITK" OFF)
endif()

set(${PRIMARY_PROJECT_NAME}_DEPENDENCIES
    SIRF
)
if (BUILD_GADGETRON)
  list(APPEND ${PRIMARY_PROJECT_NAME}_DEPENDENCIES Gadgetron)
  set(Armadillo_REQUIRED_VERSION 4.600)
endif()

if (BUILD_siemens_to_ismrmrd)
  list(APPEND ${PRIMARY_PROJECT_NAME}_DEPENDENCIES siemens_to_ismrmrd)
endif()

if (BUILD_petmr_rd_tools)
  list(APPEND ${PRIMARY_PROJECT_NAME}_DEPENDENCIES petmr_rd_tools)
endif()

ExternalProject_Include_Dependencies(${proj} DEPENDS_VAR ${PRIMARY_PROJECT_NAME}_DEPENDENCIES)

message(STATUS "")
message(STATUS "BOOST_ROOT = " ${BOOST_ROOT})
message(STATUS "ISMRMRD_DIR = " ${ISMRMRD_DIR})
message(STATUS "FFTW3_ROOT_DIR = " ${FFTW3_ROOT_DIR})
message(STATUS "STIR_DIR = " ${STIR_DIR})
message(STATUS "HDF5_ROOT = " ${HDF5_ROOT})
message(STATUS "GTEST_ROOT = " ${GTEST_ROOT})
message(STATUS "Matlab_ROOT_DIR = " ${Matlab_ROOT_DIR})
message(STATUS "PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_LIBRARIES=${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIRS=${PYTHON_INCLUDE_DIRS}")

#Need to configure main project here.
#set(proj ${PRIMARY_PROJECT_NAME})

# Make environment files
set(SIRF_SRC_DIR ${SOURCE_ROOT_DIR}/SIRF)
set(CCPPETMR_INSTALL ${SUPERBUILD_INSTALL_DIR})

## configure the environment files env_ccppetmr.sh/csh
## We create a whole bash/csh block script which does set the appropriate
## environment variables for Python and Matlab.
## in the env_ccppetmr scripts we perform a substitution of the whole block
## during the configure_file() command call below.

## Note that the MATLAB_DEST and PYTHON_DEST variables are currently set
## in External_SIRF.cmake. That's a bit confusing of course (TODO).
set(ENV_PYTHON_BASH "#####    Python not found    #####")
set(ENV_PYTHON_CSH  "#####    Python not found    #####")
if(PYTHONINTERP_FOUND)

  set (ENV_PYTHON_CSH "\
    if $?PYTHONPATH then \n\
      #setenv PYTHONPATH ${PYTHON_DEST}:$PYTHONPATH \n\
    else \n\
      #setenv PYTHONPATH ${PYTHON_DEST} \n\
      setenv SIRF_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE}")

  set (ENV_PYTHON_BASH "\
     #export PYTHONPATH=${PYTHON_DEST}:$PYTHONPATH \n\
     export SIRF_PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}")

endif()

set(ENV_MATLAB_BASH "#####     Matlab not found     #####")
set(ENV_MATLAB_CSH  "#####     Matlab not found     #####")
if (Matlab_FOUND)
  set(ENV_MATLAB_BASH "\
MATLABPATH=${MATLAB_DEST}\n\
export MATLABPATH\n\
SIRF_MATLAB_EXECUTABLE=${Matlab_MAIN_PROGRAM}\n\
export SIRF_MATLAB_EXECUTABLE")
  set(ENV_MATLAB_CSH "\
   if $?MATLABPATH then\n\
	setenv MATLABPATH ${MATLAB_DEST}:$MATLABPATH\n\
   else\n\
	setenv MATLABPATH ${MATLAB_DEST}\n\
   endif\n\
   setenv SIRF_MATLAB_EXECUTABLE ${Matlab_MAIN_PROGRAM}")
endif()

configure_file(env_ccppetmr.sh.in ${CCPPETMR_INSTALL}/bin/env_ccppetmr.sh)
configure_file(env_ccppetmr.csh.in ${CCPPETMR_INSTALL}/bin/env_ccppetmr.csh)

if(PYTHONINTERP_FOUND)
  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild/setup.py.in")
  set(SETUP_PY "${PYTHON_DEST}/setup.py")
  set(SETUP_PY_INIT "${PYTHON_DEST}/sirf/__init__.py")
  message(STATUS "setup.py: ${SETUP_PY}")

  configure_file("${SETUP_PY_IN}" "${SETUP_PY}")

  add_custom_command(OUTPUT "${SETUP_PY_INIT}"
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${PYTHON_DEST}/sirf"
    COMMAND "${CMAKE_COMMAND}" -E touch "${SETUP_PY_INIT}"
    COMMAND "${PYTHON_EXECUTABLE}" setup.py build
    DEPENDS "${SETUP_PY_IN}"
    WORKING_DIRECTORY "${PYTHON_DEST}")

  add_custom_target(pybuild_stir ALL DEPENDS ${SETUP_PY_INIT})

  install(CODE "execute_process(COMMAND\n\
    \"${PYTHON_EXECUTABLE}\" -m pip install -U -e \"${CCPPETMR_INSTALL}\")")
endif(PYTHONINTERP_FOUND)


# add tests
enable_testing()
add_test(NAME SIRF_TESTS
         COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION>
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/builds/SIRF/build/)
