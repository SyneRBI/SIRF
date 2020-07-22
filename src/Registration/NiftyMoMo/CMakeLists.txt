#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2020 University College London
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

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add source files
SET(SOURCES "BSplineTransformation.cpp" "Transformation.cpp")

# Create library
ADD_LIBRARY(NiftyMoMo ${SOURCES})

# Headers
TARGET_INCLUDE_DIRECTORIES(NiftyMoMo PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>$<INSTALL_INTERFACE:include/>"
)

# NIFTYREG
TARGET_LINK_LIBRARIES(NiftyMoMo PUBLIC ${NR_libs_full_path} ${NR_libs_full_path})
TARGET_INCLUDE_DIRECTORIES(NiftyMoMo BEFORE PUBLIC "${NIFTYREG_INCLUDE_DIRS}")

# Install
INSTALL(TARGETS NiftyMoMo DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
