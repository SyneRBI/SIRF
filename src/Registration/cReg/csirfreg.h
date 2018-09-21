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
    void* cSIRFReg_SIRFReg_open_TM(const char* filename, PTR_FLOAT ptr_TM);
    void* cSIRFReg_compose_transformations_into_single_deformation2(const void* im, const void* trans1, const void* trans2);
    void* cSIRFReg_compose_transformations_into_single_deformation3(const void* im, const void* trans1, const void* trans2, const void* trans3);
    void* cSIRFReg_compose_transformations_into_single_deformation4(const void* im, const void* trans1, const void* trans2, const void* trans3, const void* trans4);
    void* cSIRFReg_compose_transformations_into_single_deformation5(const void* im, const void* trans1, const void* trans2, const void* trans3, const void* trans4, const void* trans5);

    // NiftiImage
    void* cSIRFReg_NiftiImage_save_to_file(const void* ptr, const char* filename);
    void* cSIRFReg_NiftiImage_fill(const void* ptr, const float val);
    void* cSIRFReg_NiftiImage_deep_copy(const void* copy_ptr, const void *orig_ptr);
    void* cSIRFReg_NiftiImage_get_dimensions(const void* ptr, PTR_INT ptr_dim);
    void* cSIRFReg_NiftiImage_get_data(const void* ptr, PTR_FLOAT ptr_data);
    void* cSIRFReg_NiftiImage_maths(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type);

    // NiftiImage3D
    void* cSIRFReg_NiftiImage3D_from_PETImageData(void* ptr);
    void* cSIRFReg_NiftiImage3D_copy_data_to(const void* ptr, const void* obj);

    // NiftiImage3DTensor
    void* cSIRFReg_NiftiImage3DTensor_save_to_file_split_xyz_components(const void* ptr, const char* filename);
    void* cSIRFReg_NiftiImage3DTensor_create_from_3D_image(const void *ptr, const void* obj);
    void* cSIRFReg_NiftiImage3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr);

	// SIRFReg
    void* cSIRFReg_SIRFReg_update(void* ptr);
    void* cSIRFReg_SIRFReg_get_deformation_displacement_image(const void* ptr, const char *transform_type);

    // SIRFRegAladin methods
    void* cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(const void* ptr, const char* filename, const char* dir);
    void* cSIRFReg_SIRFReg_get_TM(const void* ptr, PTR_FLOAT ptr_TM, const char* dir);

    // SIRFRegNiftyResample
    void* cSIRFReg_SIRFRegNiftyResample_add_transformation(void* self, const void* trans, const char* type);
    void* cSIRFReg_SIRFRegNiftyResample_update(void* ptr);

    // SIRFRegImageWeightedMean
    void* cSIRFReg_SIRFRegImageWeightedMean_add_image(void* ptr, const void* obj, const float weight);
    void* cSIRFReg_SIRFRegImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight);
    void* cSIRFReg_SIRFRegImageWeightedMean_update(void* ptr);

    // SIRFRegTransformation
    void* cSIRFReg_SIRFRegTransformation_get_as_deformation_field(const void* ptr, const void* ref);
    void* cSIRFReg_SIRFRegTransformationAffine_construct_from_TM(PTR_FLOAT ptr_TM);
    void* cSIRFReg_SIRFRegTransformationDisplacement_construct_from_NiftiImage3DDisplacement(const void* ptr);
    void* cSIRFReg_SIRFRegTransformationDeformation_construct_from_NiftiImage3DDeformation(const void* ptr);

#ifndef CSIRFREG_FOR_MATLAB
}
#endif

#endif
