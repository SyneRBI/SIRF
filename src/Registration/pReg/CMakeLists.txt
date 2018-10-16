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

if(SIRFREG_BUILD_PYTHON)
  set(CMAKE_POSITION_INDEPENDENT_CODE True)

  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")

  SET_SOURCE_FILES_PROPERTIES(pysirfreg.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(pysirfreg.i PROPERTIES SWIG_FLAGS "-I${cSIRFREG_INCLUDE_DIR}")
  # find libraries and include files
  # TODO would be better to use target_include_directories

  FIND_PACKAGE(PythonLibs REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  SET(SWIG_MODULE_pysirfreg_EXTRA_DEPS "${cSIRFREG_INCLUDE_DIR}/csirfreg.h")
  SWIG_ADD_MODULE(pysirfreg python pysirfreg.i)
  SWIG_LINK_LIBRARIES(pysirfreg SIRFReg ${PYTHON_LIBRARIES})

  include_directories(${PROJECT_SOURCE_DIR}/src/common/include)

  INSTALL(TARGETS ${SWIG_MODULE_pysirfreg_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/pysirfreg.py" SIRFReg.py DESTINATION "${PYTHON_DEST}/sirf")

ADD_SUBDIRECTORY(tests)
endif(SIRFREG_BUILD_PYTHON)
