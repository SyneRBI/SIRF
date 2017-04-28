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

find_package(STIR 3.1.0)
if (STIR_FOUND)

  # need to add this explicitly to get the STIR_REGISTRIES to compile
  include_directories("${STIR_INCLUDE_DIRS}")

  set(cSTIR_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cSTIR")

  ADD_SUBDIRECTORY(cSTIR)
  ADD_SUBDIRECTORY(pSTIR)
  ADD_SUBDIRECTORY(mSTIR)

else()

  message(WARNING "STIR not found. This is probably not what you want. If you need it, set STIR_DIR.")

endif()
