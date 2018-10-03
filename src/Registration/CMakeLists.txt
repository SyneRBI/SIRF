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
OPTION(BUILD_SIRFREG_PYTHON "Build Python SWIG module" BUILD_PYTHON)
OPTION(BUILD_SIRFREG_MATLAB "Build MATLAB SWIG module" BUILD_MATLAB)

set(cSIRFREG_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cReg")

# If python or matlab desired
IF(BUILD_SIRFREG_PYTHON OR BUILD_SIRFREG_MATLAB)
	# Find Swig
	FIND_PACKAGE(SWIG REQUIRED)
	INCLUDE(${SWIG_USE_FILE})
	SET(CMAKE_SWIG_FLAGS "")
ENDIF(BUILD_SIRFREG_PYTHON OR BUILD_SIRFREG_MATLAB)

# Python
ADD_SUBDIRECTORY(pReg)
ADD_SUBDIRECTORY(mReg)
