#========================================================================
# Author: Johannes Mayer
# Copyright 2016 - 2020 University College London
# Copyright 2020 Physikalisch-Technische Bundesanstalt Berlin
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

add_library(MR_TESTS_CPP_AUXILIARY ${CMAKE_CURRENT_SOURCE_DIR}/mrtest_auxiliary_funs.cpp)
target_link_libraries(MR_TESTS_CPP_AUXILIARY PUBLIC csirf cgadgetron)

add_executable(MR_TESTS_CPLUSPLUS ${CMAKE_CURRENT_SOURCE_DIR}/mrtests.cpp)
target_link_libraries(MR_TESTS_CPLUSPLUS PUBLIC MR_TESTS_CPP_AUXILIARY csirf cgadgetron)
INSTALL(TARGETS MR_TESTS_CPLUSPLUS DESTINATION bin)

ADD_TEST(NAME MR_TESTS_CPLUSPLUS
         COMMAND MR_TESTS_CPLUSPLUS
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
