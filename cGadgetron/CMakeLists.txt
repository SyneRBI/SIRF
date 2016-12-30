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

include_directories(${FFTW3_INCLUDE_DIR})

add_library(cgt cgadgetron.cpp xgadgetron.cpp ismrmrd_fftw.cpp)

set (cGadgetron_INCLUDE_DIR "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
# copy to parent scope
#set (cGadgetron_INCLUDE_DIR "${cGadgetron_INCLUDE_DIR}" PARENT_SCOPE)

target_include_directories(cgt PUBLIC "${cGadgetron_INCLUDE_DIR}")

target_include_directories(cgt PRIVATE "${HDF5_C_INCLUDE_DIR}")
target_link_libraries(cgt iutil "${ISMRMD_LIBRARIES}" "${FFTW3_LIBRARIES}" "${HDF5_LIBRARIES}")
