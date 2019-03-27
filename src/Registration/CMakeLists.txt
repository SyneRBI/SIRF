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
# NIFTYREG
# ------------------------------------------------------------------------
FIND_PACKAGE(NIFTYREG 1.5.61)
if (NIFTYREG_FOUND)

	# Add in the NIFTYREG stuff
	INCLUDE_DIRECTORIES("${NIFTYREG_INCLUDE_DIRS}")
  	LINK_DIRECTORIES("${NIFTYREG_LIBRARY_DIRS}")
  	SET(NIFTYREG_LIBRARY_DIRS "${NIFTYREG_LIBRARY_DIRS}" PARENT_SCOPE) # need elsewhere
	
	SET(cREG_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cReg")
	INCLUDE_DIRECTORIES(${cREG_INCLUDE_DIR})
	
	ADD_SUBDIRECTORY(cReg)
	ADD_SUBDIRECTORY(pReg)
	ADD_SUBDIRECTORY(mReg)

else()
	
	message(WARNING "NIFTYREG not found. This is probably not what you want. If you need it, set NIFTYREG_DIR.")
	# Following will prevent building synergistic code 
	# (but will not change cached variable, so attempt to find NIFTYREG again on next Configure)
	set(REG_ENABLE OFF PARENT_SCOPE)

endif()