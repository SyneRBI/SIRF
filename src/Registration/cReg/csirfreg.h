/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef cSIRFREG_INTERFACE
#define cSIRFREG_INTERFACE

#ifndef CSIRFREG_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t
extern "C" {
#else
#define PTR_INT int*
#define PTR_FLOAT float*
#define PTR_DOUBLE double*
#endif

	// Common SIRFReg Object methods
	void* cSIRFReg_newObject(const char* name);
	void* cSIRFReg_objectFromFile(const char* name, const char* filename);
	void* cSIRFReg_setParameter
		(void* ptr, const char* obj, const char* name, const void* value);
	void* cSIRFReg_parameter(const void* ptr, const char* obj, const char* name);

    // SIRFRegMisc
    void* cSIRFReg_do_nifti_images_match(const void* im1, const void* im2, const float accuracy_percentage_of_max);
    void* cSIRFReg_dump_nifti_info_filename(const char* filename);
    void* cSIRFReg_dump_nifti_info_im1(const void* im1);
    void* cSIRFReg_dump_nifti_info_im2(const void* im1, const void* im2);
    void* cSIRFReg_dump_nifti_info_im3(const void* im1, const void* im2, const void* im3);
    void* cSIRFReg_dump_nifti_info_im4(const void* im1, const void* im2, const void* im3, const void* im4);
    void* cSIRFReg_dump_nifti_info_im5(const void* im1, const void* im2, const void* im3, const void* im4, const void* im5);

	// SIRFImageData
    void* cSIRFReg_SIRFImageData_from_PETImageData(void* ptr);
    void* cSIRFReg_SIRFImageData_save_to_file(const void* ptr, const char* filename);
    void* cSIRFReg_SIRFImageData_copy_data_to(const void* ptr, const void* obj);
    void* cSIRFReg_SIRFImageData_fill(const void* ptr, const float val);

	// SIRFRegImageDataDeformation
    void* cSIRFReg_SIRFImageDataDeformation_save_to_file(const void* ptr, const char* filename, const bool split_xyz);
	void* cSIRFReg_SIRFImageDataDeformation_create_from_3D_image(void* ptr, const void* obj);

	// SIRFReg
    void* cSIRFReg_SIRFReg_update(void* ptr);
    void* cSIRFReg_SIRFReg_save_image(const void* ptr, const char* filename);
    void* cSIRFReg_SIRFReg_save_deformation_displacement_image(const void* ptr, const char* filename, const char* type, const bool split_xyz);

    // SIRFRegAladin methods
    void* cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(const void* ptr, const char* filename, const char* dir);

    // SIRFRegNiftyResample
    void* cSIRFReg_SIRFRegNiftyResample_update(void* ptr);
    void* cSIRFReg_SIRFRegNiftyResample_save_resampled_image(void* ptr, const char* filename);

    // SIRFRegImageWeightedMean3D
    void* cSIRFReg_SIRFRegImageWeightedMean3D_add_image(void* ptr, const void* obj, const float weight);
    void* cSIRFReg_SIRFRegImageWeightedMean3D_add_image_filename(void* ptr, const char* filename, const float weight);
    void* cSIRFReg_SIRFRegImageWeightedMean3D_update(void* ptr);
    void* cSIRFReg_SIRFRegImageWeightedMean3D_save_image_to_file(const void* ptr, const char* filename);

    // SIRFRegImageWeightedMean4D
    void* cSIRFReg_SIRFRegImageWeightedMean4D_add_image(void* ptr, const void* obj, const float weight);
    void* cSIRFReg_SIRFRegImageWeightedMean4D_add_image_filename(void* ptr, const char* filename, const float weight);
    void* cSIRFReg_SIRFRegImageWeightedMean4D_update(void* ptr);
    void* cSIRFReg_SIRFRegImageWeightedMean4D_save_image_to_file(const void* ptr, const char* filename);

#ifndef CSIRFREG_FOR_MATLAB
}
#endif

#endif
