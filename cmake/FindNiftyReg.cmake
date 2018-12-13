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

set(NiftyReg_FOUND 0)

SET(NiftyReg_Binary_DIR CACHE PATH "NiftyReg binary directory" )

# Ask for the binary and check
IF (NOT IS_DIRECTORY ${NiftyReg_Binary_DIR})
	MESSAGE(FATAL_ERROR "Enter the binary directory of NiftyReg.")
ENDIF()

IF(EXISTS "${NiftyReg_Binary_DIR}/CMakeCache.txt")
  LOAD_CACHE("${NiftyReg_Binary_DIR}" READ_WITH_PREFIX "nifty" "NiftyReg_SOURCE_DIR" "CMAKE_INSTALL_PREFIX")
endif()

SET(NiftyReg_Source_DIR ${niftyNiftyReg_SOURCE_DIR})
SET(NiftyReg_Install_DIR ${niftyCMAKE_INSTALL_PREFIX})

# Include
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
string(REGEX REPLACE "\n$" "" NR_VERSION "${NR_VERSION}")

SET(NiftyReg_requiredLibs
    _reg_aladin
    _reg_blockMatching
    _reg_f3d
    _reg_maths
    _reg_ReadWriteImage
    _reg_resampling
    _reg_tools
    reg_nifti
    reg_png
    _reg_femTrans
    _reg_globalTrans
    _reg_localTrans
    _reg_measure
    )
 
# Loop over each library, find it and add it. 
# Do it twice to resolve linking errors depending on order
FOREACH(elem ${NiftyReg_requiredLibs} ${NiftyReg_requiredLibs})
    FIND_LIBRARY(NiftyReg_${elem} ${elem} HINTS ${NiftyReg_Install_DIR}/lib)
    MARK_AS_ADVANCED(NiftyReg_${elem})
    SET(NiftyReg_Libs ${NiftyReg_Libs} ${NiftyReg_${elem}} )
ENDFOREACH(elem ${NiftyReg_requiredLibs})

# Success!
set(NiftyReg_FOUND 1)
