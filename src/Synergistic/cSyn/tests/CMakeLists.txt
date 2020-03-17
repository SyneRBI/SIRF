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
TARGET_LINK_LIBRARIES(test_conv_img LINK_PUBLIC csirf iutilities cstir cgadgetron Reg)
INSTALL(TARGETS test_conv_img DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

set(MR_vendor_dicom_as_nifti "${CMAKE_SOURCE_DIR}/data/examples/MR/zenodo/dicom_as_nifti.nii")
set(MR_SIRF_recon "${CMAKE_SOURCE_DIR}/data/examples/MR/zenodo/SIRF_recon.h5")

ADD_TEST(NAME SYN_TEST_STIR_TO_NIFTI_CPLUSPLUS 
	COMMAND test_conv_img ${MR_vendor_dicom_as_nifti} "STIR" "Reg"
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
ADD_TEST(NAME SYN_TEST_GADGETRON_TO_NIFTI_CPLUSPLUS 
	COMMAND test_conv_img ${MR_SIRF_recon} "Gadgetron" "Reg"
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
