#========================================================================
# Author: Kris Thielemans
# Copyright 2016,2017 University College London
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

if (CMAKE_VERSION VERSION_LESS "3.1")
   set (CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
else ()
   set (CMAKE_CXX_STANDARD 11)
endif ()

# include explicit search for HDF5 for older ISMRMRD versions.
find_package(HDF5 COMPONENTS C)
find_package(ISMRMRD)
if (ISMRMRD_FOUND)
  set(CMAKE_MODULE_PATH "${ISMRMRD_DIR}")
  # We need the following line because there is a problem with the
  # ISMRMRD config cmake file which means that library location is
  # not part of the target.
  set(cGadgetron_LIBRARY_DIRS "${ISMRMRD_LIBRARY_DIRS}" PARENT_SCOPE)
  set(cGadgetron_INCLUDE_DIR "${ISMRMRD_INCLUDE_DIR}" PARENT_SCOPE)

  find_package(FFTW3 COMPONENTS single REQUIRED)

  if (SIRF_INSTALL_DEPENDENCIES AND WIN32)
    set(FFTW3_DLL_DIR ${FFTW3_INCLUDE_DIR})
    message(STATUS "Install FFTW3 libraries from ${FFTW3_DLL_DIR} ")
	file(GLOB FFTW_DLLS "${FFTW3_DLL_DIR}/*.dll")
    install( FILES ${FFTW_DLLS} DESTINATION bin )
endif()
	
    
  find_package(HDF5 1.8 COMPONENTS C REQUIRED)

  include_directories("${ISMRMRD_INCLUDE_DIR}")
  # need to add this as ISMRMRD_LIBRARIES doesn't set the path to the library
  link_directories("${ISMRMRD_LIBRARY_DIRS}")

  ADD_SUBDIRECTORY(cGadgetron)

  set(cGadgetron_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cGadgetron")

  ADD_SUBDIRECTORY(pGadgetron)
  ADD_SUBDIRECTORY(mGadgetron)

endif()
