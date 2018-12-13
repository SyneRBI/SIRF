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

set(cSIRF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

add_library(csirf csirf.cpp)
target_include_directories(csirf PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>"
  )
target_link_libraries(csirf iutilities)

# TODO target
include_directories(include)
add_library(common GeometricalInfo.cpp)
#INSTALL(common)

if (BUILD_PYTHON)

  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")
  SET_SOURCE_FILES_PROPERTIES(pysirf.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(pysirf.i PROPERTIES SWIG_FLAGS "-I${cSIRF_INCLUDE_DIR}")
  FIND_PACKAGE(PythonLibs REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  set(SWIG_MODULE_pysirf_EXTRA_DEPS "csirf.h")
  SWIG_ADD_MODULE(pysirf python pysirf.i)
  SWIG_LINK_LIBRARIES(pysirf csirf iutilities ${PYTHON_LIBRARIES})
  INSTALL(TARGETS ${SWIG_MODULE_pysirf_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/pysirf.py" DESTINATION "${PYTHON_DEST}/sirf")
  file(GLOB PythonFiles "${CMAKE_CURRENT_LIST_DIR}/*.py")
  INSTALL(FILES ${PythonFiles} DESTINATION "${PYTHON_DEST}/sirf")

endif (BUILD_PYTHON)

if (BUILD_MATLAB)

  add_library(mig mig.cpp)

  INSTALL(DIRECTORY +mUtilities DESTINATION "${MATLAB_DEST}")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

endif (BUILD_MATLAB)
