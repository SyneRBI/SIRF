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

if(BUILD_PYTHON)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13") 
    # policy introduced in CMake 3.13
    cmake_policy(SET CMP0078 OLD)
  endif()
  set(CMAKE_POSITION_INDEPENDENT_CODE True)

  # Modify STIR.py depending on presence of NiftyPET
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/STIR.py.in" "${CMAKE_CURRENT_BINARY_DIR}/STIR.py" @ONLY)

  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")

  SET(cSTIR_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cSTIR/include")

  SET_SOURCE_FILES_PROPERTIES(pystir.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(pystir.i PROPERTIES SWIG_FLAGS "-I${cSTIR_INCLUDE_DIR}")
  # find libraries and include files
  # TODO would be better to use target_include_directories

  FIND_PACKAGE(PythonLibs REQUIRED)
  # Tell CMake it needs to re-run SWIG when .h file changes
  set(SWIG_MODULE_pystir_EXTRA_DEPS "${cSTIR_INCLUDE_DIR}/sirf/STIR/cstir.h")
  SWIG_ADD_LIBRARY(pystir LANGUAGE python TYPE MODULE SOURCES pystir.i  ${STIR_REGISTRIES})
  TARGET_INCLUDE_DIRECTORIES(${SWIG_MODULE_pystir_REAL_NAME} PUBLIC ${PYTHON_INCLUDE_PATH})
  SWIG_LINK_LIBRARIES(pystir csirf cstir iutilities ${STIR_LIBRARIES} ${PYTHON_LIBRARIES})

  INSTALL(TARGETS ${SWIG_MODULE_pystir_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES 
    "${CMAKE_CURRENT_BINARY_DIR}/pystir.py"
    "${CMAKE_CURRENT_BINARY_DIR}/STIR.py"
    "${CMAKE_CURRENT_SOURCE_DIR}/STIR_params.py"
    DESTINATION "${PYTHON_DEST}/sirf")

ADD_SUBDIRECTORY(tests)
endif(BUILD_PYTHON)
