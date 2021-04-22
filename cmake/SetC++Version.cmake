#========================================================================
# Author: Kris Thielemans
# Copyright 2020-2021 University College London
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

# A macro to set the C++ version
# If CMAKE_CXX_STANDARD is already set to a more recent version, keep that.
macro(UseCXX VERSION)
  if (NOT DEFINED CMAKE_CXX_STANDARD)
    set (CMAKE_CXX_STANDARD ${VERSION})
  else()
    if ((CMAKE_CXX_STANDARD EQUAL 98) OR (CMAKE_CXX_STANDARD LESS VERSION))
      message(FATAL_ERROR "CXX_STANDARD needs to be at least ${VERSION}")
    endif()
  endif ()
endmacro(UseCXX)
