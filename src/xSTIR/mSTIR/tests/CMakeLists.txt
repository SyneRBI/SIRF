#========================================================================
# Author: Kris Thielemans
# Copyright 2017 University College London
#
# This file is part of the CCP PETMR Synergistic Image Reconstruction Framework (SIRF) SuperBuild.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#=========================================================================

# next commands work (on Windows) but would need further work to get the log file handled properly
# Luckily CMake's FindMatlab provides a function for this. We'll use that below.

#set(test_cmd "pwd;ver;try; test_all; catch err, fprintf(2,\'ERROR: %s\\\\n\', err.message),exit(1), end, exit(0)")
#message(STATUS "test_cmd ${test_cmd}")
#add_test(NAME PET_TESTS_MATLAB
#         COMMAND matlab -nodesktop -wait -logfile "${CMAKE_CURRENT_BINARY_DIR}/test_all.log" -r "${test_cmd}"
#         WORKING_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR})

matlab_add_unit_test(NAME PET_TESTS_MATLAB
    UNITTEST_FILE test_all.m
    ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}
    NO_UNITTEST_FRAMEWORK)

# Only test the NiftyPET python wrapper if STIR has been built with NiftyPET
if (STIR_WITH_NIFTYPET_PROJECTOR)
	matlab_add_unit_test(NAME PET_TESTS_MATLAB_NIFTYPET
	    UNITTEST_FILE test_niftypet.m
	    ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}
	    NO_UNITTEST_FRAMEWORK)
endif()