#========================================================================
# Author: Richard Brown
# Copyright 2018 University College London
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

# Create an executable
ADD_EXECUTABLE (SYNERGISTIC_TEST_CPLUSPLUS test_cSynergistic.cpp ${STIR_REGISTRIES})
SET_TARGET_PROPERTIES (SYNERGISTIC_TEST_CPLUSPLUS
    PROPERTIES FOLDER ${CMAKE_INSTALL_PREFIX}/bin
    INSTALL_NAME_DIR "Synergistic")
TARGET_LINK_LIBRARIES(SYNERGISTIC_TEST_CPLUSPLUS PUBLIC csirf iutilities cstir cgadgetron SIRFReg)

# Make into test
ADD_TEST(NAME SYNERGISTIC_TEST_CPLUSPLUS COMMAND SYNERGISTIC_TEST_CPLUSPLUS ${CMAKE_SOURCE_DIR} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


######################################################################
# For now, just return since there are no matlab and python files
######################################################################
return()

# If python
if (BUILD_PYTHON)
	# Make into test
	ADD_TEST(NAME SYNERGISTIC_TEST_PYTHON
	COMMAND ${PYTHON_EXECUTABLE} -m nose ${CMAKE_CURRENT_SOURCE_DIR}/test_pSynergistic.py
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

# If matlab
if (BUILD_MATLAB)
	MATLAB_ADD_UNIT_TEST(NAME SYNERGISTIC_TEST_MATLAB
	UNITTEST_FILE   test_mSynergistic.m
	ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}
	NO_UNITTEST_FRAMEWORK)
endif()