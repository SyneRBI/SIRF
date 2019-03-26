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
FIND_PACKAGE(NiftyReg 1.5.61)
if (NiftyReg_FOUND)

	get_cmake_property(_variableNames VARIABLES)
	list (SORT _variableNames)
	foreach (_variableName ${_variableNames})
	    message(STATUS "${_variableName}=${${_variableName}}")
	endforeach()

	# Add in the NiftyReg stuff
	INCLUDE_DIRECTORIES("${NiftyReg_INCLUDE_DIRS}")
  	LINK_DIRECTORIES("${NiftyReg_LIBRARY_DIRS}")
  	SET(NiftyReg_LIBRARY_DIRS "${NiftyReg_LIBRARY_DIRS}" PARENT_SCOPE) # need elsewhere
	
	SET(cREG_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cReg")
	INCLUDE_DIRECTORIES(${cREG_INCLUDE_DIR})
	
	ADD_SUBDIRECTORY(cReg)
	ADD_SUBDIRECTORY(pReg)
	ADD_SUBDIRECTORY(mReg)

else()

  message(WARNING "NiftyReg not found. This is probably not what you want. If you need it, set NiftyReg_DIR.")

endif()