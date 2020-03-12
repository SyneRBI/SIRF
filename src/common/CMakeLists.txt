#========================================================================
# Author: Evgueni Ovtchinnikov
# Copyright 2017 University College London
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
target_link_libraries(csirf PUBLIC iutilities)
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
  SWIG_ADD_LIBRARY(pysirf LANGUAGE python TYPE MODULE SOURCES pysirf.i)
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
  add_library(msirf SHARED msirf.c)
  target_include_directories(msirf PUBLIC ${Matlab_INCLUDE_DIRS})
  SET_TARGET_PROPERTIES(msirf PROPERTIES SUFFIX ".${MATLAB_MEX_EXT}" PREFIX "${MATLAB_PREFIX}") 
  target_link_libraries(msirf csirf iutilities ${Matlab_LIBRARIES})

  INSTALL(TARGETS msirf DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES msirf.h DESTINATION "${MATLAB_DEST}")
  INSTALL(DIRECTORY +sirf DESTINATION "${MATLAB_DEST}")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

  ADD_SUBDIRECTORY(gmi)

endif (BUILD_MATLAB)
