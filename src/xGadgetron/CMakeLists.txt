#========================================================================
# Author: Kris Thielemans, Richard Brown
# Copyright 2016-2019 University College London
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

if (SIRF_INSTALL_DEPENDENCIES AND WIN32)
  set(FFTW3_DLL_DIR ${FFTW3_INCLUDE_DIR})
  message(STATUS "Install FFTW3 libraries from ${FFTW3_DLL_DIR} ")
  file(GLOB FFTW_DLLS "${FFTW3_DLL_DIR}/*.dll")
  install( FILES ${FFTW_DLLS} DESTINATION bin )
endif()

ADD_SUBDIRECTORY(cGadgetron)
ADD_SUBDIRECTORY(pGadgetron)
ADD_SUBDIRECTORY(mGadgetron)