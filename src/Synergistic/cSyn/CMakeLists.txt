#========================================================================
# Author: Richard Brown
# Author: Evgueni Ovtchinnikov
# Copyright 2018 - 2020 University College London
# Copyright 2020 UKRI-STFC
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

SET(CMAKE_POSITION_INDEPENDENT_CODE True)

#add_library(csyn csyn.cpp utilities.cpp)
add_library(csyn utilities.cpp)
target_include_directories(csyn PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>"
  )
target_include_directories(csyn PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>"
  )
target_link_libraries(csyn PUBLIC cgadgetron cstir Reg)
target_compile_definitions(csyn PRIVATE WIN32_LEAN_AND_MEAN)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
  COMPONENT Development
  DESTINATION include)
install(TARGETS csyn EXPORT SIRFTargets
  COMPONENT Development
  DESTINATION lib)

ADD_SUBDIRECTORY(tests)
