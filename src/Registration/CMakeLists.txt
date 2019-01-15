#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2019 University College London
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

# ------------------------------------------------------------------------
# NiftyReg
# ------------------------------------------------------------------------
FIND_PACKAGE(NiftyReg)
if (NiftyReg_FOUND)
	if(NOT ("${NR_VERSION}" MATCHES "1.5.59" OR "${NR_VERSION}" MATCHES "1.5.61"))
		MESSAGE(STATUS "\n\nSIRF's registration was developed with NiftyReg 1.5.59 and 1.5.61, and cannot be guaranteed for other version numbers.\n"
	    "Your compiled version of NiftyReg is version ${NR_VERSION}.\n"
	    "If there is a version mismatch, the keys in the parser may need to be altered.\n")
	endif()
else()
	MESSAGE(FATAL_ERROR "Registration has been enabled, but there are no registration packages available.")
endif()

# Include a preprocessor to compile code depenging on version number.
add_definitions(-DNR_VERSION="${NR_VERSION}")

# ------------------------------------------------------------------------
# Add c++, python, matlab components (if desired)
# ------------------------------------------------------------------------
SET(cREG_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cReg")
INCLUDE_DIRECTORIES(${cREG_INCLUDE_DIR})
ADD_SUBDIRECTORY(cReg)
ADD_SUBDIRECTORY(pReg)
ADD_SUBDIRECTORY(mReg)
