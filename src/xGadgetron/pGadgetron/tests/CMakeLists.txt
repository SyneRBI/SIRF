#========================================================================
# Copyright 2017 - 2021 Science Technology Facilities Council
# Copyright 2017 - 2022 University College London
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
add_test(NAME MR_TESTS_PYTHON
  COMMAND ${Python_EXECUTABLE} -m pytest --cov=sirf.Gadgetron --cov-config=.coveragerc-Gadgetron src/xGadgetron/pGadgetron/
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
add_test(NAME MR_DEMOS_PYTHON
  COMMAND ${Python_EXECUTABLE} run_all.py --non-interactive
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/Python/MR")
add_test(NAME MR_GADGETRON_DEMOS_PYTHON
  COMMAND ${Python_EXECUTABLE} run_all.py --non-interactive
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/examples/Python/MR/Gadgetron")
