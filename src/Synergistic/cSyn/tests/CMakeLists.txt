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
#  limitations under the License.
#
#=========================================================================

SET(CMAKE_POSITION_INDEPENDENT_CODE True)

ADD_EXECUTABLE(test_conv_img test_conv_img.cpp ${STIR_REGISTRIES})
TARGET_LINK_LIBRARIES(test_conv_img LINK_PUBLIC csirf iutilities cstir cgadgetron Reg utilities)
INSTALL(TARGETS test_conv_img DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
