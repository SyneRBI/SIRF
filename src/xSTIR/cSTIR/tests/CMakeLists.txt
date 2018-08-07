#========================================================================
# Author: Kris Thielemans
# Copyright 2016 University College London
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
#ADD_EXECUTABLE (test4 ${CMAKE_CURRENT_SOURCE_DIR}/test4.cpp ${STIR_REGISTRIES})
#SET_TARGET_PROPERTIES (test4
#    PROPERTIES FOLDER ${CMAKE_INSTALL_PREFIX}/bin
#    INSTALL_NAME_DIR "${PROJECT_NAME}")
#TARGET_LINK_LIBRARIES(test4 PUBLIC cstir )

  include_directories(${PROJECT_SOURCE_DIR}/src/common/include)
  add_executable(test4 ${CMAKE_CURRENT_SOURCE_DIR}/test4.cpp ${STIR_REGISTRIES})
  target_link_libraries(test4 cstir ${STIR_LIBRARIES})
  INSTALL(TARGETS test4 DESTINATION bin)

ADD_TEST(NAME PET_TEST_CPLUSPLUS COMMAND test4 WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
