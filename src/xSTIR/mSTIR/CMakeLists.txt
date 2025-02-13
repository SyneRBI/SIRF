#========================================================================
# Author: Kris Thielemans, Richard Brown
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

  set(CMAKE_POSITION_INDEPENDENT_CODE True)

  set(mstir_c ${CMAKE_CURRENT_BINARY_DIR}/mstir.c)
  set(mstir_h ${CMAKE_CURRENT_BINARY_DIR}/mstir.h)

  add_executable(gmi_xstir gmi/gmi_xstir.cpp)
  target_link_libraries(gmi_xstir csirf mig)
  # Generate m-files on the fly
  add_custom_command(
    OUTPUT ${mstir_c} ${mstir_h}
    COMMAND gmi_xstir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    DEPENDS ${CMAKE_SOURCE_DIR}/src/xSTIR/cSTIR/include/sirf/STIR/cstir.h
    COMMENT "Generating mstir.c/.h"
    VERBATIM
  )

  add_library(mstir SHARED ${mstir_c} printer.cpp  ${STIR_REGISTRIES})
  target_include_directories(mstir PUBLIC ${Matlab_INCLUDE_DIRS})
  # fix output name and link flags
  SET_TARGET_PROPERTIES(mstir PROPERTIES
        SUFFIX ".${MATLAB_MEX_EXT}" PREFIX "${MATLAB_PREFIX}") 
  target_link_libraries(mstir csirf cstir iutilities ${STIR_LIBRARIES} ${Matlab_LIBRARIES})
  INSTALL(TARGETS mstir DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES ${mstir_h} DESTINATION "${MATLAB_DEST}")
  INSTALL(DIRECTORY +sirf DESTINATION "${MATLAB_DEST}")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

  ADD_SUBDIRECTORY(tests)
endif(BUILD_MATLAB)
