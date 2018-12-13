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

# NiftyReg
find_package(NiftyReg)
if (NOT NiftyReg_FOUND)
	message(FATAL_ERROR "NiftyReg not found.")
endif()

MESSAGE(STATUS "\n\nSIRFReg was developed with NiftyReg 1.5.59, and cannot be guaranteed for other version numbers.\n"
    "Your compiled version of NiftyReg is version ${NR_VERSION}.\n"
    "If there is a version mismatch, the keys in the parser may need to be altered.\n")

# Compile the C++
ADD_SUBDIRECTORY(cReg)
INCLUDE_DIRECTORIES(cReg)

set(cSIRFREG_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cReg")

ADD_SUBDIRECTORY(pReg)
ADD_SUBDIRECTORY(mReg)
