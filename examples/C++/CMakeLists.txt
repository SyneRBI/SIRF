#========================================================================
# Author: Kris Thielemans
# Copyright 2021 University College London
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

# cmake file for using SIRF. See the SIRF User's Guide and http://www.cmake.org.

cmake_minimum_required(VERSION 3.9.0)

PROJECT(SIRFTest)

# SIRF exports a SIRFConfig.cmake etc. This defines targets for every
# SIRF library, e.g. sirf::common, sirf::Reg etc
# Let's find SIRF (when running CMake, you might have to set SIRF_DIR
# to let CMake find the exported SIRF).
find_package(SIRF 3.1.0)

# Add an executable that uses a part of SIRF

add_executable(example_using_sirf example_using_sirf.cpp)
# Currently need to explicitly pass STIR_LIBRARIES due to a problem in STIR
target_link_libraries(example_using_sirf sirf::cstir ${STIR_LIBRARIES})
