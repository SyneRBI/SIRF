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

find_package(Boost COMPONENTS system thread filesystem REQUIRED)

link_directories(${ISMRMRD_LIBRARY_DIRS})

if(BUILD_PYTHON)

  FIND_PACKAGE(SWIG REQUIRED)
  INCLUDE("${SWIG_USE_FILE}")

  SET_SOURCE_FILES_PROPERTIES(pygadgetron.i PROPERTIES CPLUSPLUS ON)  
  SET_SOURCE_FILES_PROPERTIES(pygadgetron.i PROPERTIES SWIG_FLAGS "-I${cGadgetron_INCLUDE_DIR}")
  # find libraries and include files
  # TODO would be better to use target_include_directories

  FIND_PACKAGE(PythonLibs REQUIRED)
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
  #FIND_PACKAGE(Numpy REQUIRED)
  #INCLUDE_DIRECTORIES(${NUMPY_INCLUDE_DIRS})
  SWIG_ADD_MODULE(pygadgetron python pygadgetron.i)
  SWIG_LINK_LIBRARIES(pygadgetron cgt iutil ismrmrd fftw3 ${ISMRMRD_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}) #fftw3f

  INSTALL(TARGETS ${SWIG_MODULE_pygadgetron_REAL_NAME} DESTINATION ${PYTHON_DEST})
  INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/pygadgetron.py  DESTINATION ${PYTHON_DEST})

endif(BUILD_PYTHON)
