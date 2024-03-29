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

	OPTION(DISABLE_SPM "disable SPM" OFF)
	if (NOT DISABLE_SPM)
		find_package(SPM QUIET)
		if(SPM_FOUND)
			MESSAGE(STATUS "SPM found")
			set(SPM_BOOL_STR "1")
		else()
			MESSAGE(STATUS "SPM not found")
			set(SPM_BOOL_STR "0")
		endif()
	else()
		SET(SPM_FOUND OFF)
		set(SPM_BOOL_STR "0")
	endif()
else()
	MESSAGE(STATUS "Matlab not enabled, so SPM cannot be enabled")
	set(SPM_BOOL_STR "0")
endif()
set(SPM_BOOL_STR ${SPM_BOOL_STR} PARENT_SCOPE)

# Check for VTK
OPTION(DISABLE_VTK "disable VTK" OFF)
if(NOT DISABLE_VTK)
	find_package(VTK QUIET)
	MESSAGE(STATUS "VTK found: ${VTK_FOUND}")
        set(SIRF_BUILT_WITH_VTK TRUE PARENT_SCOPE)
        set(VTK_VERSION ${VTK_VERSION} PARENT_SCOPE)
else()
	set(VTK_FOUND OFF)
endif()

add_subdirectory(NiftyMoMo)
add_subdirectory(cReg)
if (BUILD_PYTHON)
  add_subdirectory(pReg)
endif()
add_subdirectory(mReg)
