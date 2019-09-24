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
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")

  include_directories(${PROJECT_SOURCE_DIR}/src/common/include)
  include_directories(${PROJECT_SOURCE_DIR}/src/xGadgetron/cGadgetron/include)

  SET_SOURCE_FILES_PROPERTIES(pygadgetron.i PROPERTIES CPLUSPLUS ON)  
  SET_SOURCE_FILES_PROPERTIES(pygadgetron.i PROPERTIES SWIG_FLAGS "-I${cGadgetron_INCLUDE_DIR}")
  # find libraries and include files
  # TODO would be better to use target_include_directories

  FIND_PACKAGE(PythonLibs REQUIRED)
  INCLUDE_DIRECTORIES("${PYTHON_INCLUDE_PATH}")
  #FIND_PACKAGE(Numpy REQUIRED)
  #INCLUDE_DIRECTORIES("${NUMPY_INCLUDE_DIRS}")
  # Tell CMake it needs to re-run SWIG when .h file changes
  set(SWIG_MODULE_pygadgetron_EXTRA_DEPS "${cGadgetron_INCLUDE_DIR}/cgadgetron.h")
  SWIG_ADD_LIBRARY(pygadgetron LANGUAGE python TYPE MODULE SOURCES pygadgetron.i)
  SWIG_LINK_LIBRARIES(pygadgetron cgadgetron iutilities "${PYTHON_LIBRARIES}")

  INSTALL(TARGETS ${SWIG_MODULE_pygadgetron_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/pygadgetron.py" Gadgetron.py DESTINATION "${PYTHON_DEST}/sirf")

  SET(LOWERCASE_PACKAGE_NAME "pygadgetron")
  CONFIGURE_FILE("${CMAKE_SOURCE_DIR}/src/common/parameters.py.in" "${CMAKE_CURRENT_BINARY_DIR}/params/Gadgetron.py")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/params/Gadgetron.py" DESTINATION "${PYTHON_DEST}/sirf/params")

ADD_SUBDIRECTORY(tests)
endif(BUILD_PYTHON)

