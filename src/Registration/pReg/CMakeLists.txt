#========================================================================
# Author: Richard Brown
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

<<<<<<< HEAD
# Set p_or_m to Python for the configure file
SET(p_or_m p)

=======
>>>>>>> ed941c275db9c90a8093c65c9490027d67adab0c
# Find python
FIND_PACKAGE( PythonLibs REQUIRED )
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

<<<<<<< HEAD
CONFIGURE_FILE( ../swigDoxy/SWIG_SIRFReg.in ${CMAKE_CURRENT_BINARY_DIR}/p${PROJECT_NAME}.i )

# Swig to convert c++ to python
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/p${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_LIBRARY(p${PROJECT_NAME} LANGUAGE python SOURCES ${CMAKE_CURRENT_BINARY_DIR}/p${PROJECT_NAME}.i)
SWIG_LINK_LIBRARIES(p${PROJECT_NAME} ${PYTHON_LIBRARIES} ${PROJECT_NAME})

# Install
INSTALL(TARGETS  						 _p${PROJECT_NAME}	  DESTINATION ${PYTHON_DEST})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/p${PROJECT_NAME}.py DESTINATION ${PYTHON_DEST})
=======
# Swig to convert c++ to python
SET_SOURCE_FILES_PROPERTIES(p${PROJECT_NAME}.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_LIBRARY(p${PROJECT_NAME} LANGUAGE python SOURCES p${PROJECT_NAME}.i)
SWIG_LINK_LIBRARIES(p${PROJECT_NAME} ${PYTHON_LIBRARIES} ${PROJECT_NAME})

# Install
INSTALL(TARGETS  						 _p${PROJECT_NAME}	  DESTINATION ${CMAKE_INSTALL_PREFIX}/python)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/p${PROJECT_NAME}.py DESTINATION ${CMAKE_INSTALL_PREFIX}/python)
>>>>>>> ed941c275db9c90a8093c65c9490027d67adab0c

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(test)
ENDIF(BUILD_TESTING)