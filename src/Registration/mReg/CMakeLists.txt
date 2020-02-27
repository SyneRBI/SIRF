#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2020 University College London
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

  SET(SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/Transformation.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/Registration.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/setParameter.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftyResample.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/parameter.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/Quaternion.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftiImageData3DTensor.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftyAladinSym.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftyF3dSym.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftyRegistration.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftiImageData3DDeformation.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftiImageData3DDisplacement.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/ImageWeightedMean.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftiImageData.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/NiftiImageData3D.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/AffineTransformation.m"
    "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/ImageData.m")

  if (SPM_FOUND)
    set(SOURCES "${SOURCES}" "${CMAKE_CURRENT_SOURCE_DIR}/+sirf/+Reg/SPMRegistration.m")
  endif()

  add_library(mreg SHARED mreg.c printer.cpp)
  target_include_directories(mreg PUBLIC ${Matlab_INCLUDE_DIRS})
  # fix output name and link flags
  SET_TARGET_PROPERTIES(mreg PROPERTIES
        SUFFIX ".${MATLAB_MEX_EXT}" PREFIX "${MATLAB_PREFIX}") 
  target_link_libraries(mreg Reg iutilities csirf cstir ${Matlab_LIBRARIES})
  INSTALL(TARGETS mreg DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES mreg.h DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES ${SOURCES} DESTINATION "${MATLAB_DEST}/+sirf/+Reg")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

  ADD_SUBDIRECTORY(gmi)
  ADD_SUBDIRECTORY(tests)
endif(BUILD_MATLAB)
