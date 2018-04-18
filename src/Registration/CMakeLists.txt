#========================================================================
# Author: Richard Brown
# Copyright 2016, 2017 University College London
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

PROJECT(SIRFReg)

# NiftyReg
SET(NiftyReg_Source_DIR  CACHE PATH "NiftyReg source directory" )
SET(NiftyReg_Binary_DIR  CACHE PATH "NiftyReg binary directory" )
SET(NiftyReg_Install_DIR CACHE PATH "NiftyReg install directory")

# To be able to build with niftyreg, you need to append the cmake_cxx_flags
# to suppress errors from converting float to double
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

# Include STIR
FIND_PACKAGE(STIR 3.1.0)
IF(STIR_FOUND)
    INCLUDE_DIRECTORIES("${STIR_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/xSTIR/cSTIR)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/iUtilities)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/common/include)
ENDIF(STIR_FOUND)

IF ((NOT IS_DIRECTORY ${NiftyReg_Source_DIR}) AND (NOT IS_DIRECTORY ${NiftyReg_Binary_DIR}) AND (NOT IS_DIRECTORY ${NiftyReg_Install_DIR}))
	MESSAGE(FATAL_ERROR "Enter the source, binary and install directories of NiftyReg.")
ENDIF()

INCLUDE_DIRECTORIES(${NiftyReg_Source_DIR}/reg-io/nrrd/NrrdIO) # ugly, but required for nrrdIO.h
INCLUDE_DIRECTORIES(${NiftyReg_Binary_DIR})                    # ugly, but required for NrrdConfigure.h
INCLUDE_DIRECTORIES(${NiftyReg_Install_DIR}/include)

# Get the NiftyReg versions from its CMakeLists.txt
# This for older versions of NiftyReg
FILE(READ ${NiftyReg_Source_DIR}/CMakeLists.txt NiftyReg_CMakeLists)
STRING(REGEX MATCH "NiftyReg_VERSION_MAJOR ([0-9]*)" _ ${NiftyReg_CMakeLists} )
SET(NiftyReg_VERSION_MAJOR ${CMAKE_MATCH_1})
STRING(REGEX MATCH "NiftyReg_VERSION_MINOR ([0-9]*)" _ ${NiftyReg_CMakeLists} )
SET(NiftyReg_VERSION_MINOR ${CMAKE_MATCH_1})
STRING(REGEX MATCH "NiftyReg_VERSION_PATCH ([0-9]*)" _ ${NiftyReg_CMakeLists} )
SET(NiftyReg_VERSION_PATCH ${CMAKE_MATCH_1})
# This for more recent versions of NiftyReg
IF (NOT DEFINED NiftyReg_VERSION_MAJOR)
	STRING(REGEX MATCH "NR_VERSION_MAJOR ([0-9]*)" _ ${NiftyReg_CMakeLists} )
	SET(NiftyReg_VERSION_MAJOR ${CMAKE_MATCH_1})
	STRING(REGEX MATCH "NR_VERSION_MINOR ([0-9]*)" _ ${NiftyReg_CMakeLists} )
	SET(NiftyReg_VERSION_MINOR ${CMAKE_MATCH_1})
	FILE(READ ${NiftyReg_Source_DIR}/niftyreg_build_version.txt NiftyReg_VERSION_PATCH)
ENDIF()
SET(NR_VERSION "${NiftyReg_VERSION_MAJOR}.${NiftyReg_VERSION_MINOR}.${NiftyReg_VERSION_PATCH}")
MESSAGE(STATUS "\n\nSIRFReg was developed with NiftyReg 1.5.0, and cannot be guaranteed for other version numbers.")
MESSAGE(STATUS "Your compiled version of NiftyReg is version ${NR_VERSION}.")
MESSAGE(STATUS "If there is a version mismatch, the keys in the parser may need to be altered.\n")
IF (NR_VERSION LESS "1.4")
	add_definitions(-DNIFTYREG_VER_1_3)
ELSEIF(NR_VERSION LESS "1.6")
	add_definitions(-DNIFTYREG_VER_1_5)
ENDIF()

# Compile the C++
ADD_SUBDIRECTORY(cReg)
INCLUDE_DIRECTORIES(cReg)

# Python or MATLAB
OPTION(BUILD_PYTHON "Build Python SWIG module" OFF)
OPTION(BUILD_MATLAB "Build MATLAB SWIG module" OFF)

# If python or matlab desired
IF(BUILD_PYTHON OR BUILD_MATLAB)
	# Find Swig
	FIND_PACKAGE(SWIG REQUIRED)
	INCLUDE(${SWIG_USE_FILE})
	SET(CMAKE_SWIG_FLAGS "")

	# Wrap the doxygen documentation in matlab and python
	OPTION(SIRFREG_SWIG_DOXY "Wrap the SIRFReg doxygen documentation so that it is available in Matlab and Python." ON)
	IF(SIRFREG_SWIG_DOXY)
		ADD_SUBDIRECTORY(swigDoxy)
		SET(SWIG_SIRFREG_DOCUMENTATION_PATH "%include ${CMAKE_CURRENT_BINARY_DIR}/swigDoxy/SWIG_SIRFREG_DOCUMENTATION.i")
	ENDIF(SIRFREG_SWIG_DOXY)
ENDIF(BUILD_PYTHON OR BUILD_MATLAB)

# Python
IF(BUILD_PYTHON)
    ADD_SUBDIRECTORY(pReg)
ENDIF()

# MATLAB
IF(BUILD_MATLAB)
    ADD_SUBDIRECTORY(mReg)
ENDIF()
