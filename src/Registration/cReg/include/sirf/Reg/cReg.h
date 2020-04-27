/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2017 - 2020 University College London

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

#ifndef cREG_INTERFACE
#define cREG_INTERFACE

#ifndef CREG_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t

extern "C" {
#else
#define PTR_INT int*
#define PTR_FLOAT float*
#define PTR_DOUBLE double*
#endif

    // Common Reg Object methods
    void* cReg_newObject(const char* name);
    void* cReg_objectFromFile(const char* name, const char* filename);
    void* setParameter(void* ptr, const char* obj, const char* name, const void* value);
    void* parameter(const void* ptr, const char* obj, const char* name);

    // NiftiImageData
    void* cReg_NiftiImageData_print_headers(const void *handle_vector_ptr);
    void* cReg_NiftiImageData_write(const void* ptr, const char* filename, const int datatype);
    void* cReg_NiftiImageData_fill(const void* ptr, const float val);
    void* cReg_NiftiImageData_fill_arr(const void* ptr, PTR_FLOAT val);
    void* cReg_NiftiImageData_deep_copy(const void* copy_ptr, const void *orig_ptr);
    void* cReg_NiftiImageData_get_dimensions(const void* ptr, PTR_INT ptr_dim);
    void* cReg_NiftiImageData_get_voxel_sizes(const void* ptr, PTR_FLOAT ptr_out);
    void* cReg_NiftiImageData_as_array(const void* ptr, PTR_FLOAT ptr_data);
    void* cReg_NiftiImageData_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type);
    void* cReg_NiftiImageData_maths_num(const void *res_ptr, const void* im1_ptr, const float val, const int maths_type);
    void* cReg_NiftiImageData_equal(const void* im1_ptr, const void* im2_ptr);
    void* cReg_NiftiImageData_norm(const void* im1_ptr, const void* im2_ptr);
    void* cReg_NiftiImageData_get_original_datatype(const void* im_ptr);
    void* cReg_NiftiImageData_crop(const void* im_ptr, PTR_INT min_index_ptr, PTR_INT max_index_ptr);
    void* cReg_NiftiImageData_pad(const void* im_ptr, PTR_INT min_index_ptr, PTR_INT max_index_ptr, const float val);
    void* cReg_NiftiImageData_set_voxel_spacing(const void* im_ptr, const float x, const float y, const float z, const int interpolation_order);
    void* cReg_NiftiImageData_normalise_zero_and_one(const void* im_ptr);
    void* cReg_NiftiImageData_standardise(const void* im_ptr);
    void* cReg_NiftiImageData_get_inner_product(const void* im1_ptr, const void* im2_ptr);
    void* cReg_NiftiImageData_from_SIRFImageData(void* ptr);
    void* cReg_NiftiImageData_from_complex_ImageData_real_component(void* in_ptr);
    void* cReg_NiftiImageData_from_complex_ImageData_imag_component(void* in_ptr);
    void* cReg_NiftiImageData_are_equal_to_given_accuracy(void* im1_ptr, void* im2_ptr, const float accuracy);

    // NiftiImageData3D

    // NiftiImageData3DTensor
    void* cReg_NiftiImageData3DTensor_write_split_xyz_components(const void* ptr, const char* filename, const int datatype);
    void* cReg_NiftiImageData3DTensor_create_from_3D_image(const void *ptr, const void* obj);
    void* cReg_NiftiImageData3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr);
    void* cReg_NiftiImageData3DTensor_flip_component(const void *ptr, const int dim);

    // NiftiImageData3DDeformation
    void* cReg_NiftiImageData3DDeformation_compose_single_deformation(const void* im, const char* types, const void* trans_vector_ptr);
    void* cReg_NiftiImageData3DDeformation_create_from_disp(const void* disp_ptr);
    void* cReg_NiftiImageData3DDeformation_get_inverse(const void* def_ptr, const void* floating_ptr);

    // NiftiImageData3DDisplacement
    void* cReg_NiftiImageData3DDisplacement_create_from_def(const void* def_ptr);

    // Registration
    void* cReg_Registration_process(void* ptr);
    void* cReg_Registration_get_deformation_displacement_image(const void* ptr, const char *transform_type, const int idx);
    void* cReg_Registration_add_floating(const void* ptr, const void *im_ptr);
    void* cReg_Registration_clear_floatings(const void* ptr);
    void* cReg_Registration_get_output(const void* ptr,const int idx);
    void* cReg_Registration_set_reference_image_filename(const void* ptr, const char* filename);
    void* cReg_Registration_set_floating_image_filename(const void* ptr, const char* filename);
    void* cReg_Registration_add_floating_image_filename(const void* ptr, const char* filename);

    // NiftyReg-based registration
    void* cReg_NiftyRegistration_set_parameter(const void* ptr, const char* par, const char* arg1, const char* arg2);
    void* cReg_NiftyRegistration_print_all_wrapped_methods(const char* name);

    // Aladin methods
    void* cReg_NiftyAladin_get_TM(const void* ptr, const char* dir);

    // SPM methods
    void* cReg_SPMRegistration_get_TM(const void* ptr, const char* dir, const int idx);

    // NiftyResample
    void* cReg_NiftyResample_add_transformation(void* self, const void* trans, const char* type);
    void* cReg_NiftyResample_clear_transformations(void* self);
    void* cReg_NiftyResample_process(void* ptr);
    void* cReg_NiftyResample_forward(const void *output_ptr, const void * const input_ptr, const void *resampler_ptr);
    void* cReg_NiftyResample_adjoint(const void *output_ptr, const void * const input_ptr, const void *resampler_ptr);

    // ImageWeightedMean
    void* cReg_ImageWeightedMean_add_image(void* ptr, const void* obj, const float weight);
    void* cReg_ImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight);
    void* cReg_ImageWeightedMean_process(void* ptr);

    // Transformation
    void* cReg_Transformation_get_as_deformation_field(const void* ptr, const char* name, const void* ref);

    // AffineTransformation
    void* cReg_AffineTransformation_construct_from_TM(PTR_FLOAT ptr_TM);
    void* cReg_AffineTransformation_construct_from_trans_and_quaternion(PTR_FLOAT trans_ptr, const void* quat_ptr);
    void* cReg_AffineTransformation_construct_from_trans_and_euler(PTR_FLOAT trans_ptr, PTR_FLOAT euler_ptr);
    void* cReg_AffineTransformation_deep_copy(const void* ptr);
    void* cReg_AffineTransformation_write(const void* ptr, const char* filename);
    void* cReg_AffineTransformation_as_array(const void* ptr, PTR_FLOAT ptr_TM);
    void* cReg_AffineTransformation_get_identity();
    void* cReg_AffineTransformation_get_inverse(const void* ptr);
    void* cReg_AffineTransformation_get_Euler_angles(const void* ptr, PTR_FLOAT Euler);
    void* cReg_AffineTransformation_get_quaternion(const void* ptr);
    void* cReg_AffineTransformation_mul(const void* mat1_ptr, const void* mat2_ptr);
    void* cReg_AffineTransformation_equal(const void* mat1_ptr, const void* mat2_ptr);
    void* cReg_AffineTransformation_get_average(const void *handle_vector_ptr);

    // Quaternion
    void* cReg_Quaternion_construct_from_array(PTR_FLOAT arr);
    void* cReg_Quaternion_construct_from_AffineTransformation(const void* ptr);
    void* cReg_Quaternion_get_average(const void *handle_vector_ptr);
    void* cReg_Quaternion_as_array(const void* ptr, PTR_FLOAT arr);

#ifndef CREG_FOR_MATLAB
}
#endif

#endif
