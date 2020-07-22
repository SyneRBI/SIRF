#========================================================================
# Author: Evgueni Ovtchinnikov
# Copyright 2020 University College London
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

set(CMAKE_POSITION_INDEPENDENT_CODE True)

set(cSIRF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_library(csirf csirf.cpp ImageData.cpp GeometricalInfo.cpp)
target_include_directories(csirf PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>"
  )
target_include_directories(csirf PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>"
  )
if ((NOT DISABLE_STIR) AND (NOT DISABLE_Gadgetron) AND (NOT DISABLE_Registration) AND "${STIR_BUILT_WITH_ITK}")
  MESSAGE(STATUS "Registration, ISMRMRD and STIR (with ITK) have been built.")
  target_link_libraries(csirf PUBLIC iutilities csyn)
else()
  MESSAGE(STATUS "Either Registration or ISMRMRD or STIR (with ITK) have not been built.")
  add_library(Syn syn_utilities.cpp)
  target_include_directories(Syn PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>"
  )
  target_include_directories(Syn PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>"
  )
  target_link_libraries(Syn PUBLIC cgadgetron cstir Reg)
  INSTALL(TARGETS Syn DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  target_link_libraries(csirf PUBLIC iutilities Syn)
endif()
INSTALL(TARGETS csirf DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

if (BUILD_PYTHON)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13") 
    # policy introduced in CMake 3.13
    cmake_policy(SET CMP0078 OLD)
  endif()
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")
  SET_SOURCE_FILES_PROPERTIES(pysirf.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(pysirf.i PROPERTIES SWIG_FLAGS "-I${cSIRF_INCLUDE_DIR}")
  FIND_PACKAGE(PythonLibs REQUIRED)
  set(SWIG_MODULE_pysirf_EXTRA_DEPS "${cSIRF_INCLUDE_DIR}/sirf/common/csirf.h")
  SWIG_ADD_LIBRARY(pysirf LANGUAGE python TYPE MODULE SOURCES pysirf.i ${STIR_REGISTRIES})
  TARGET_INCLUDE_DIRECTORIES(${SWIG_MODULE_pysirf_REAL_NAME} PUBLIC ${PYTHON_INCLUDE_PATH})
  SWIG_LINK_LIBRARIES(pysirf csirf iutilities ${PYTHON_LIBRARIES})
  INSTALL(TARGETS ${SWIG_MODULE_pysirf_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/pysirf.py" DESTINATION "${PYTHON_DEST}/sirf")
  file(GLOB PythonFiles "${CMAKE_CURRENT_LIST_DIR}/*.py")
  INSTALL(FILES ${PythonFiles} DESTINATION "${PYTHON_DEST}/sirf")

endif (BUILD_PYTHON)

if (BUILD_MATLAB)

  add_library(mig mig.cpp)
  target_include_directories(mig PUBLIC ${Matlab_INCLUDE_DIRS})

  set(msirf_c ${CMAKE_CURRENT_BINARY_DIR}/msirf.c)
  set(msirf_h ${CMAKE_CURRENT_BINARY_DIR}/msirf.h)

  add_executable(gmi_csirf gmi/gmi_csirf.cpp)
  target_link_libraries(gmi_csirf mig csirf)
  # Generate m-files on the fly
  add_custom_command(
    OUTPUT ${msirf_c} ${msirf_h}
    COMMAND gmi_csirf ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    DEPENDS ${CMAKE_SOURCE_DIR}/src/common/include/sirf/common/csirf.h
    COMMENT "Generating msirf.c/.h"
    VERBATIM
  )

  add_library(msirf SHARED ${msirf_c})
  target_include_directories(msirf PUBLIC ${Matlab_INCLUDE_DIRS})
  SET_TARGET_PROPERTIES(msirf PROPERTIES SUFFIX ".${MATLAB_MEX_EXT}" PREFIX "${MATLAB_PREFIX}") 
  target_link_libraries(msirf csirf iutilities ${Matlab_LIBRARIES})

  INSTALL(TARGETS msirf DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES ${msirf_h} DESTINATION "${MATLAB_DEST}")
  INSTALL(DIRECTORY +sirf DESTINATION "${MATLAB_DEST}")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

endif (BUILD_MATLAB)
