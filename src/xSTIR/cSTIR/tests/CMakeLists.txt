#========================================================================
# Author: Kris Thielemans
# Copyright 2016, 2018, 2021 University College London
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
  add_executable(cstir_tests test1.cpp test2.cpp main.cpp ${STIR_REGISTRIES})
  target_link_libraries(cstir_tests csirf cstir ${STIR_LIBRARIES})
  INSTALL(TARGETS cstir_tests DESTINATION bin)

  add_executable(cstir_test4  test4.cpp run_test4.cpp ${STIR_REGISTRIES})
  target_link_libraries(cstir_test4 csirf cstir ${STIR_LIBRARIES})
  INSTALL(TARGETS cstir_test4 DESTINATION bin)

  add_executable(cstir_test6  test6.cpp run_test6.cpp ${STIR_REGISTRIES})
  target_link_libraries(cstir_test6 csirf cstir ${STIR_LIBRARIES})
  INSTALL(TARGETS cstir_test6 DESTINATION bin)

ADD_TEST(NAME PET_TESTS_CPLUSPLUS_1 COMMAND cstir_tests WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

ADD_TEST(NAME PET_TESTS_CPLUSPLUS_4 COMMAND cstir_test4 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#ADD_TEST(NAME PET_TESTS_CPLUSPLUS_6 COMMAND cstir_test6 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

