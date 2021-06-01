#========================================================================
# Author: Johannes Mayer
# Copyright 2016 - 2020 University College London
# Copyright 2021 Physikalisch-Technische Bundesanstalt Berlin
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

set( SIMULATION_TEST_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_dynamicsimulation.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_tissueparameters.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_contrastgenerator.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_phantom_input.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_dynamics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_noisegenerator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_volume_orientator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_dynsim_deformer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_auxiliary_testing_functions.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests_auxiliary_input_output.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/auxiliary_testing_functions.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/app_mracquisitiondata.cpp
)

add_library(SIMULATION_TESTS_LIB ${SIMULATION_TEST_FILES})
target_link_libraries(SIMULATION_TESTS_LIB PUBLIC csirf cgadgetron cdynamicsimulation MR_TESTS_CPP_AUXILIARY)
target_include_directories(SIMULATION_TESTS_LIB PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")


add_executable(SIMULATION_TESTS_CPLUSPLUS ${CMAKE_CURRENT_SOURCE_DIR}/all_simulation_tests.cpp ${STIR_REGISTRIES})
target_link_libraries(SIMULATION_TESTS_CPLUSPLUS PUBLIC SIMULATION_TESTS_LIB ${STIR_LIBRARIES})
INSTALL(TARGETS SIMULATION_TESTS_CPLUSPLUS DESTINATION bin)

ADD_TEST(NAME SIMULATION_TESTS_CPLUSPLUS
         COMMAND SIMULATION_TESTS_CPLUSPLUS
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
