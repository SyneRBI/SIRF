#========================================================================
# Author: Edoardo Pasca
# Copyright 2017 University College London
# Copyright 2017 Science Technology Facilities Council
# Author: Casper da Costa-Luis
# Copyright 2017
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

# Only test the NiftyPET python wrapper if STIR has been built with NiftyPET
if (STIR_WITH_NIFTYPET_PROJECTOR)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/tests_niftypet.py.in ${CMAKE_CURRENT_SOURCE_DIR}/tests_niftypet.py COPYONLY)
endif()

add_test(NAME PET_TESTS_PYTHON
  COMMAND ${PYTHON_EXECUTABLE} -m nose --with-coverage --cover-package=sirf.STIR src/xSTIR/pSTIR/
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
