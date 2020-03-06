#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2020 University College London
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

# Have we got SPM? (only works with Matlab)
if (BUILD_MATLAB)
	# if on linux, engine is in funny place
	if("${CMAKE_SYSTEM}" MATCHES "Linux")
		list(APPEND Matlab_LIBRARIES 
			"${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so" 
			"${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so")
	endif()

	find_package(SPM QUIET)
	if(SPM_FOUND)
		MESSAGE(STATUS "SPM found")
		set(SPM_BOOL_STR "1")
	else()
		MESSAGE(STATUS "SPM not found")
		set(SPM_BOOL_STR "0")
	endif()
else()
	MESSAGE(STATUS "Matlab not enabled, so SPM cannot be enabled")
	set(SPM_BOOL_STR "0")
endif()

ADD_SUBDIRECTORY(NiftyMoMo)
ADD_SUBDIRECTORY(cReg)
ADD_SUBDIRECTORY(pReg)
ADD_SUBDIRECTORY(mReg)