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

add_library(cgt cgadgetron.cpp xgadgetron.cpp ismrmrd_fftw.cpp)

set (cGadgetron_INCLUDE_DIR "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>$<INSTALL_INTERFACE:include>")
# copy to parent scope
#set (cGadgetron_INCLUDE_DIR "${cGadgetron_INCLUDE_DIR}" PARENT_SCOPE)

target_include_directories(cgt PUBLIC "${cGadgetron_INCLUDE_DIR}")

target_include_directories(cgt PRIVATE "${HDF5_C_INCLUDE_DIR}")
target_link_libraries(cgt iutil "${ISMRMD_LIBRARIES}" "${FFTW3_LIBRARIES}" "${HDF5_LIBRARIES}")
