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

if(BUILD_PYTHON)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13") 
    # policy introduced in CMake 3.13
    cmake_policy(SET CMP0078 OLD)
  endif()
  set(CMAKE_POSITION_INDEPENDENT_CODE True)

  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")

  # Modify pReg.py depending on presence of SPM
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Reg.py.in" "${CMAKE_CURRENT_BINARY_DIR}/Reg.py")

  SET(cREG_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cReg/include")

  SET_SOURCE_FILES_PROPERTIES(pyreg.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(pyreg.i PROPERTIES SWIG_FLAGS "-I${cREG_INCLUDE_DIR}")
  # find libraries and include files
  # TODO would be better to use target_include_directories

  # Tell CMake it needs to re-run SWIG when .h file changes
  set(SWIG_MODULE_pyreg_EXTRA_DEPS "${cREG_INCLUDE_DIR}/sirf/Reg/cReg.h")

  FIND_PACKAGE(PythonLibs REQUIRED)
  SWIG_ADD_LIBRARY(pyreg LANGUAGE python TYPE MODULE SOURCES pyreg.i)
  TARGET_INCLUDE_DIRECTORIES(${SWIG_MODULE_pyreg_REAL_NAME} PUBLIC ${PYTHON_INCLUDE_PATH})
  SWIG_LINK_LIBRARIES(pyreg PUBLIC Reg ${PYTHON_LIBRARIES})

  INSTALL(TARGETS ${SWIG_MODULE_pyreg_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/pyreg.py"
    "${CMAKE_CURRENT_BINARY_DIR}/Reg.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/Reg_params.py"
    DESTINATION "${PYTHON_DEST}/sirf")

ADD_SUBDIRECTORY(tests)
endif(BUILD_PYTHON)
