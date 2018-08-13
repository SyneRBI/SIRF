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

# To be able to build with niftyreg, you need to append the cmake_cxx_flags
# to suppress errors from converting float to double
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

# SIRF stuff
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/common/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/iUtilities)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/common/include/SIRF/common)
# STIR stuff
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/xSTIR/cSTIR)
find_package(STIR)
INCLUDE_DIRECTORIES("${STIR_INCLUDE_DIRS}")
# Gadgetron stuff
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/xGadgetron/cGadgetron)
find_package(ISMRMRD)
INCLUDE_DIRECTORIES(${ISMRMRD_INCLUDE_DIRS})

# NiftyReg
include("${CMAKE_CURRENT_SOURCE_DIR}/NiftyReg.cmake")

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
#    ADD_SUBDIRECTORY(pReg)
ENDIF()

# MATLAB
IF(BUILD_MATLAB)
#    ADD_SUBDIRECTORY(mReg)
ENDIF()
