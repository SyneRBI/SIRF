#========================================================================
# Author: Richard Brown
# Copyright 2018 - 2019 University College London
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





# add_library( NiftyMoMoTransformation 
#                BSplineTransformation.h 
#                BSplineTransformation.cpp 
#                Transformation.h
# 			   Transformation.cpp )

# Add source files
SET(SOURCES "BSplineTransformation.cpp" "Transformation.cpp")

# Create library
ADD_LIBRARY(NiftyMoMo ${SOURCES})

# Headers
TARGET_INCLUDE_DIRECTORIES(NiftyMoMo PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>$<INSTALL_INTERFACE:include/>"
)

# If niftyreg was bulit with OpenMP
if (NIFTYREG_BUILT_WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  if (OpenMP_CXX_FOUND)
    TARGET_LINK_LIBRARIES(NiftyMoMo PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()

# Install
INSTALL(TARGETS NiftyMoMo DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
INSTALL (
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
