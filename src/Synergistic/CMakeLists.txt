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

# Need this because there's a problem in ISMRMRDConfig.cmake. 
# This has to come before making any executables
LINK_DIRECTORIES(${cGadgetron_LIBRARY_DIRS})
INCLUDE_DIRECTORIES(${cGadgetron_INCLUDE_DIR})

# NIFTYREG libraries
LINK_DIRECTORIES("${NIFTYREG_LIBRARY_DIRS}")

# Executables
SET(SYN_executables
    sirf_registration
    sirf_do_images_match
    sirf_resample
    )

FOREACH(elem ${SYN_executables})
    ADD_EXECUTABLE(${elem} ${elem}.cpp)
    TARGET_LINK_LIBRARIES(${elem} LINK_PUBLIC csirf iutilities cstir cgadgetron Reg)
    INSTALL(TARGETS ${elem} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
ENDFOREACH(elem ${SYN_executables})

# Tests
ADD_SUBDIRECTORY(tests)