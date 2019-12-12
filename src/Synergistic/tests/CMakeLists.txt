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

function(get_zenodo_file path url md5)
	message(STATUS "Downloading (if not already done)...
		Path:     ${path}
		URL:      ${url}
		MD5:      ${md5}")
	file(DOWNLOAD "${url}" "${path}" SHOW_PROGRESS EXPECTED_MD5 "${md5}")
endfunction()

OPTION(DOWNLOAD_ZENODO_TEST_DATA "download zenodo data for tests" ON)
if (DOWNLOAD_ZENODO_TEST_DATA)
	set(MR_SIRF_recon "${CMAKE_SOURCE_DIR}/data/examples/MR/zenodo/SIRF_recon.h5")
	get_zenodo_file(
		"${MR_SIRF_recon}"
		"https://zenodo.org/record/3571228/files/SIRF_recon.h5?download=1"
		"0d69892502d3a5cd220716a717cc2c26")
	set(MR_vendor_dicom_as_nifti "${CMAKE_SOURCE_DIR}/data/examples/MR/zenodo/dicom_as_nifti.nii")
	get_zenodo_file(
		"${MR_vendor_dicom_as_nifti}"
		"https://zenodo.org/record/3571228/files/dicom_as_nifti.nii?download=1"
		"c5f5b519e5c032560ead912713479b53")
endif()

# Create an executable
ADD_EXECUTABLE (SYN_TEST_CPLUSPLUS test_cSynergistic.cpp ${STIR_REGISTRIES})
SET_TARGET_PROPERTIES (SYN_TEST_CPLUSPLUS
    PROPERTIES FOLDER ${CMAKE_INSTALL_PREFIX}/bin
    INSTALL_NAME_DIR "Synergistic")
TARGET_LINK_LIBRARIES(SYN_TEST_CPLUSPLUS PUBLIC csirf iutilities cstir cgadgetron Reg)

if (DOWNLOAD_ZENODO_TEST_DATA)
	set(test_data "${MR_vendor_dicom_as_nifti}" "${MR_SIRF_recon}")
else()
	set(test_data "${CMAKE_SOURCE_DIR}/data/examples/Registration/test2.nii.gz")
endif()

# Make into test
ADD_TEST(NAME SYN_TEST_CPLUSPLUS 
	COMMAND SYN_TEST_CPLUSPLUS ${test_data}
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# If python
if (BUILD_PYTHON)
	# Make into test
	ADD_TEST(NAME SYN_TEST_PYTHON
	COMMAND ${PYTHON_EXECUTABLE} -m nose ${CMAKE_CURRENT_SOURCE_DIR}/test_pSynergistic.py
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

# If matlab
if (BUILD_MATLAB)
	MATLAB_ADD_UNIT_TEST(NAME SYN_TEST_MATLAB
	UNITTEST_FILE   test_mSynergistic.m
	ADDITIONAL_PATH ${CMAKE_CURRENT_SOURCE_DIR}
	NO_UNITTEST_FRAMEWORK
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()