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

########################################################################################
# test SIRFReg
########################################################################################
# Create an executable
ADD_EXECUTABLE (SIRFREG_TEST_CPLUSPLUS test_cSIRFReg.cpp ${STIR_REGISTRIES})
SET_TARGET_PROPERTIES (SIRFREG_TEST_CPLUSPLUS
    PROPERTIES FOLDER ${CMAKE_INSTALL_PREFIX}/bin
    INSTALL_NAME_DIR "SIRFReg")
TARGET_LINK_LIBRARIES(SIRFREG_TEST_CPLUSPLUS PUBLIC SIRFReg)

# Make into test
ADD_TEST(NAME SIRFREG_TEST_CPLUSPLUS COMMAND SIRFREG_TEST_CPLUSPLUS WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#[[
########################################################################################
# test f3d
########################################################################################
ADD_EXECUTABLE (SIRFREG_TEST_F3D f3d_mouse.cpp)
SET_TARGET_PROPERTIES (SIRFREG_TEST_F3D
    PROPERTIES FOLDER ${CMAKE_INSTALL_PREFIX}/bin
    INSTALL_NAME_DIR "SIRFReg")
TARGET_LINK_LIBRARIES(SIRFREG_TEST_F3D PUBLIC SIRFReg )

# Make into test
ADD_TEST(NAME SIRFREG_TEST_F3D COMMAND SIRFREG_TEST_F3D WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

########################################################################################
# Delete below after debugging
########################################################################################
ADD_EXECUTABLE (SIRFREG_TEST_ALADIN aladin_longitudinal.cpp)
SET_TARGET_PROPERTIES (SIRFREG_TEST_ALADIN
    PROPERTIES FOLDER ${CMAKE_INSTALL_PREFIX}/bin
    INSTALL_NAME_DIR "SIRFReg")
TARGET_LINK_LIBRARIES(SIRFREG_TEST_ALADIN PUBLIC SIRFReg )

# Make into test
ADD_TEST(NAME SIRFREG_TEST_ALADIN COMMAND SIRFREG_TEST_ALADIN WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})]]
