#========================================================================
# Copyright 2017 - 2021 University College London
# Copyright 2017 - 2021 Science Technology Facilities Council
#
# This file is part of the SyneRBI Synergistic Image Reconstruction Framework (SIRF) SuperBuild.
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
set (exclude_tests "--ignore-files=\"test_algebra\"" "--exclude=TestDataContainerAlgebra")
# If not NiftyPET projector, dont test that
if (NiftyPET_BOOL_STR EQUAL "0")
	set(exclude_tests "--ignore-files=\"test_algebra|tests_NiftyPET\"" "--exclude=TestDataContainerAlgebra")
endif()

add_test(NAME PET_TESTS_PYTHON
  COMMAND ${PYTHON_EXECUTABLE} -m nose --with-coverage --cover-package=sirf.STIR src/xSTIR/pSTIR ${exclude_tests}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
#add_test(NAME PET_DEMOS_PYTHON
#  COMMAND ${PYTHON_EXECUTABLE} run_all.py --non-interactive
#  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/Python/PET")
add_test(NAME PET_PYTHON_ALGEBRA
  COMMAND ${PYTHON_EXECUTABLE} -m unittest test_algebra
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})