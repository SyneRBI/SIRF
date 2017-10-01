#========================================================================
# Author: Edoardo Pasca
# Copyright 2017 University College London
# Copyright 2017 Science Technology Facilities Council
#
# This file is part of the CCP PETMR Synergistic Image Reconstruction Framework (SIRF) SuperBuild.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#         http://www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#=========================================================================

add_test(NAME MR_FULLY_SAMPLED
         COMMAND ${PYTHON_EXECUTABLE} fully_sampled.py WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/xGadgetron/pGadgetron/tests/ )

add_test(NAME MR_UNDER_SAMPLED
         COMMAND ${PYTHON_EXECUTABLE} undersampled.py WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/xGadgetron/pGadgetron/tests/ )

