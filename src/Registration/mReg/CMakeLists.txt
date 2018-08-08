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

# Set p_or_m to MATLAB for the configure file
SET(p_or_m m)

# Find matlab
FIND_PACKAGE 		 (Matlab COMPONENTS MAIN_PROGRAM)
INCLUDE_DIRECTORIES  (${Matlab_INCLUDE_DIRS})

# libmx isn't set by find_package(matlab), so set it ourselves
# Be careful! There is a libmx.dylib in /usr/lib (on OSX, at least),
# so make sure you link to the matlab version (they are not the same!)
FIND_LIBRARY		 (libmx mx HINTS ${Matlab_BINARIES_DIR})

CONFIGURE_FILE( ../SWIG_SIRFReg.in ${CMAKE_CURRENT_BINARY_DIR}/mSIRFReg.i )

# Swig to convert c++ to matlab
SET_PROPERTY	     (SOURCE ${CMAKE_CURRENT_BINARY_DIR}/mSIRFReg.i PROPERTY CPLUSPLUS ON)
SWIG_ADD_LIBRARY     (mSIRFReg LANGUAGE matlab SOURCES ${CMAKE_CURRENT_BINARY_DIR}/mSIRFReg.i ${STIR_REGISTRIES})
SWIG_LINK_LIBRARIES  (mSIRFReg SIRFReg #[[SIRF_common]] ${Matlab_LIBRARIES} ${libmx})
SET_TARGET_PROPERTIES(mSIRFReg PROPERTIES SUFFIX "MEX.${Matlab_MEX_EXTENSION}" )

# Install
INSTALL(TARGETS mSIRFReg DESTINATION ${MATLAB_DEST})
INSTALL(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/+mSIRFReg DESTINATION ${MATLAB_DEST})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigGet.m DESTINATION ${MATLAB_DEST})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigRef.m DESTINATION ${MATLAB_DEST})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/SwigMem.m DESTINATION ${MATLAB_DEST})

ADD_SUBDIRECTORY(tests)