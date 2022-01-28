#========================================================================
# Author: Kris Thielemans
# Copyright 2017 - 2020 University College London
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

if(BUILD_MATLAB)
  #commented out as done at top-level
  #find_package(Boost COMPONENTS system thread filesystem REQUIRED)

  set(CMAKE_POSITION_INDEPENDENT_CODE True)

  set(mgadgetron_c ${CMAKE_CURRENT_BINARY_DIR}/mgadgetron.c)
  set(mgadgetron_h ${CMAKE_CURRENT_BINARY_DIR}/mgadgetron.h)

  add_executable(gmi_xgadgetron gmi/gmi_xgadgetron.cpp)
  target_link_libraries(gmi_xgadgetron mig csirf)
  # Generate m-files on the fly
  add_custom_command(
    OUTPUT ${mgadgetron_c} ${mgadgetron_h}
    COMMAND gmi_xgadgetron ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    DEPENDS ${CMAKE_SOURCE_DIR}/src/xGadgetron/cGadgetron/include/sirf/Gadgetron/cgadgetron.h
    COMMENT "Generating mgadgetron.c/.h"
    VERBATIM
  )

  add_library(mgadgetron SHARED ${mgadgetron_c})
  target_include_directories(mgadgetron PUBLIC ${Matlab_INCLUDE_DIRS})
  # fix output name and link flags
  SET_TARGET_PROPERTIES(mgadgetron PROPERTIES
        SUFFIX ".${MATLAB_MEX_EXT}" PREFIX "${MATLAB_PREFIX}") 
  target_link_libraries(mgadgetron csirf cgadgetron iutilities ${Matlab_LIBRARIES} )

  INSTALL(TARGETS mgadgetron DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES ${mgadgetron_h} DESTINATION "${MATLAB_DEST}")
  INSTALL(DIRECTORY +sirf DESTINATION "${MATLAB_DEST}")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

  ADD_SUBDIRECTORY(tests)
endif(BUILD_MATLAB)
