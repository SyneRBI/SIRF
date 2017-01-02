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

set(CMAKE_POSITION_INDEPENDENT_CODE True)

if (CMAKE_VERSION VERSION_LESS "3.1")
   set (CMAKE_CXX_FLAGS "--std=c++11 ${CMAKE_CXX_FLAGS}")
else ()
   set (CMAKE_CXX_STANDARD 11)
endif ()
  
find_package(ISMRMRD REQUIRED)
set(CMAKE_MODULE_PATH ${ISMRMRD_DIR})
find_package(FFTW3 COMPONENTS single double REQUIRED)

find_package(HDF5 1.8 COMPONENTS C REQUIRED)

include_directories(${ISMRMRD_INCLUDE_DIR})
include_directories(${FFTW3_INCLUDE_DIR})

ADD_SUBDIRECTORY(cGadgetron)
ADD_SUBDIRECTORY(pGadgetron)
