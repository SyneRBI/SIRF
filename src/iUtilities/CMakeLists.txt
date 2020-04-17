#========================================================================
# Author: Kris Thielemans
# Copyright 2016, 2017 University College London
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

add_library(iutilities iutilities.cpp)
target_include_directories(iutilities PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>"
  )
target_include_directories(iutilities PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>$<INSTALL_INTERFACE:include>"
  )

target_include_directories(iutilities PUBLIC ${Boost_INCLUDE_DIRS})

if(BUILD_PYTHON)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.13") 
    # policy introduced in CMake 3.13
    cmake_policy(SET CMP0078 OLD)
  endif()
  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")

  SET(ciUtilities_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")

  SET_SOURCE_FILES_PROPERTIES(pyiutilities.i PROPERTIES CPLUSPLUS ON)
  SET_SOURCE_FILES_PROPERTIES(pyiutilities.i PROPERTIES SWIG_FLAGS "-I${ciUtilities_INCLUDE_DIR}")
  # find libraries and include files
  # TODO would be better to use target_include_directories

  FIND_PACKAGE(PythonLibs REQUIRED)
  set(SWIG_MODULE_pyiutilities_EXTRA_DEPS 
    "${ciUtilities_INCLUDE_DIR}/sirf/iUtilities/DataHandle.h" 
    "${ciUtilities_INCLUDE_DIR}/sirf/iUtilities/iutilities.h")
  SWIG_ADD_LIBRARY(pyiutilities LANGUAGE python TYPE MODULE SOURCES pyiutilities.i)
  TARGET_INCLUDE_DIRECTORIES(${SWIG_MODULE_pyiutilities_REAL_NAME} PUBLIC ${PYTHON_INCLUDE_PATH})
  SWIG_LINK_LIBRARIES(pyiutilities iutilities ${PYTHON_LIBRARIES})

  INSTALL(TARGETS ${SWIG_MODULE_pyiutilities_REAL_NAME} DESTINATION "${PYTHON_DEST}/sirf")
  INSTALL(FILES "${CMAKE_CURRENT_BINARY_DIR}/pyiutilities.py"  DESTINATION "${PYTHON_DEST}/sirf")

endif(BUILD_PYTHON)

if(BUILD_MATLAB)

  include_directories(${Matlab_INCLUDE_DIRS})

  add_library(miutilities SHARED miutilities.c)
  # fix output name and link flags
  SET_TARGET_PROPERTIES(miutilities PROPERTIES
        SUFFIX ".${MATLAB_MEX_EXT}" PREFIX "${MATLAB_PREFIX}") 
  target_link_libraries(miutilities iutilities ${Matlab_LIBRARIES} )

  INSTALL(TARGETS miutilities DESTINATION "${MATLAB_DEST}")
  INSTALL(FILES miutilities.h DESTINATION "${MATLAB_DEST}")

  ADD_SUBDIRECTORY(gmi)
endif(BUILD_MATLAB)
