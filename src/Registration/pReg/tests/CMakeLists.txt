#========================================================================
# Author: Richard Brown
# Copyright 2016-2020 University College London
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
if (SPM_BOOL_STR EQUAL "0")
  set (exclude_tests "--ignore-files=\"test_algebra\"" "--exclude=TestDataContainerAlgebra")
else()
  set (exclude_tests "--ignore-files=\"test_algebra\"" "--exclude=TestDataContainerAlgebra")
endif()
# Make into test
ADD_TEST(NAME REG_TESTS_PYTHON
  COMMAND ${Python_EXECUTABLE} -m nose --with-coverage --cover-package=sirf.Reg --cover-erase -d -v ${CMAKE_SOURCE_DIR}/src/Registration/pReg ${exclude_tests})

add_test(NAME REG_PYTHON_ALGEBRA
  COMMAND ${Python_EXECUTABLE} -m unittest test_algebra
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})