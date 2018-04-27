#========================================================================
# Author: Evgueni Ovtchinnikov
# Copyright 2017 University College London
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

if (BUILD_PYTHON)

  file(GLOB PythonFiles "${CMAKE_CURRENT_LIST_DIR}/*.py")
  INSTALL(FILES ${PythonFiles} DESTINATION "${PYTHON_DEST}/sirf")

endif (BUILD_PYTHON)

if (BUILD_MATLAB)

  add_library(mig mig.cpp)

  INSTALL(DIRECTORY +mUtilities DESTINATION "${MATLAB_DEST}")
  file(GLOB MatlabFiles "${CMAKE_CURRENT_LIST_DIR}/*.m")
  INSTALL(FILES ${MatlabFiles} DESTINATION "${MATLAB_DEST}")

endif (BUILD_MATLAB)

