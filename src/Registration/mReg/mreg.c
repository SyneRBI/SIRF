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
#define CREG_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#include <mex.h>
#include "matrix.h"
#include "sirf/Reg/cReg.h"

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
EXPORTED_FUNCTION  void* mReg_newObject(const char* name) {
	return cReg_newObject(name);
}
EXPORTED_FUNCTION     void* mReg_objectFromFile(const char* name, const char* filename) {
	return cReg_objectFromFile(name, filename);
}
EXPORTED_FUNCTION     void* mSetParameter(void* ptr, const char* obj, const char* name, const void* value) {
	return setParameter(ptr, obj, name, value);
}
EXPORTED_FUNCTION     void* mParameter(const void* ptr, const char* obj, const char* name) {
	return parameter(ptr, obj, name);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_print_headers(const void *handle_vector_ptr) {
	return cReg_NiftiImageData_print_headers(handle_vector_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_write(const void* ptr, const char* filename, const int datatype) {
	return cReg_NiftiImageData_write(ptr, filename, datatype);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_fill(const void* ptr, const float val) {
	return cReg_NiftiImageData_fill(ptr, val);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_fill_arr(const void* ptr, PTR_FLOAT val) {
	return cReg_NiftiImageData_fill_arr(ptr, val);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_deep_copy(const void* copy_ptr, const void *orig_ptr) {
	return cReg_NiftiImageData_deep_copy(copy_ptr, orig_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_dimensions(const void* ptr, PTR_INT ptr_dim) {
	return cReg_NiftiImageData_get_dimensions(ptr, ptr_dim);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_voxel_sizes(const void* ptr, PTR_FLOAT ptr_out) {
	return cReg_NiftiImageData_get_voxel_sizes(ptr, ptr_out);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_as_array(const void* ptr, PTR_FLOAT ptr_data) {
	return cReg_NiftiImageData_as_array(ptr, ptr_data);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type) {
	return cReg_NiftiImageData_maths_im(res_ptr, im1_ptr, im2_ptr, maths_type);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_maths_num(const void *res_ptr, const void* im1_ptr, const float val, const int maths_type) {
	return cReg_NiftiImageData_maths_num(res_ptr, im1_ptr, val, maths_type);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_equal(const void* im1_ptr, const void* im2_ptr) {
	return cReg_NiftiImageData_equal(im1_ptr, im2_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_norm(const void* im1_ptr, const void* im2_ptr) {
	return cReg_NiftiImageData_norm(im1_ptr, im2_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_original_datatype(const void* im_ptr) {
	return cReg_NiftiImageData_get_original_datatype(im_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_crop(const void* im_ptr, PTR_INT min_index_ptr, PTR_INT max_index_ptr) {
	return cReg_NiftiImageData_crop(im_ptr, min_index_ptr, max_index_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_pad(const void* im_ptr, PTR_INT min_index_ptr, PTR_INT max_index_ptr, const float val) {
	return cReg_NiftiImageData_pad(im_ptr, min_index_ptr, max_index_ptr, val);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_set_voxel_spacing(const void* im_ptr, const float x, const float y, const float z, const int interpolation_order) {
	return cReg_NiftiImageData_set_voxel_spacing(im_ptr, x, y, z, interpolation_order);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_normalise_zero_and_one(const void* im_ptr) {
	return cReg_NiftiImageData_normalise_zero_and_one(im_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_standardise(const void* im_ptr) {
	return cReg_NiftiImageData_standardise(im_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_inner_product(const void* im1_ptr, const void* im2_ptr) {
	return cReg_NiftiImageData_get_inner_product(im1_ptr, im2_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_from_SIRFImageData(void* ptr) {
	return cReg_NiftiImageData_from_SIRFImageData(ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_from_complex_ImageData_real_component(void* in_ptr) {
	return cReg_NiftiImageData_from_complex_ImageData_real_component(in_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_from_complex_ImageData_imag_component(void* in_ptr) {
	return cReg_NiftiImageData_from_complex_ImageData_imag_component(in_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData_are_equal_to_given_accuracy(void* im1_ptr, void* im2_ptr, const float accuracy) {
	return cReg_NiftiImageData_are_equal_to_given_accuracy(im1_ptr, im2_ptr, accuracy);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_write_split_xyz_components(const void* ptr, const char* filename, const int datatype) {
	return cReg_NiftiImageData3DTensor_write_split_xyz_components(ptr, filename, datatype);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_create_from_3D_image(const void *ptr, const void* obj) {
	return cReg_NiftiImageData3DTensor_create_from_3D_image(ptr, obj);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr) {
	return cReg_NiftiImageData3DTensor_construct_from_3_components(obj, x_ptr, y_ptr, z_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_flip_component(const void *ptr, const int dim) {
	return cReg_NiftiImageData3DTensor_flip_component(ptr, dim);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDeformation_compose_single_deformation(const void* im, const char* types, const void* trans_vector_ptr) {
	return cReg_NiftiImageData3DDeformation_compose_single_deformation(im, types, trans_vector_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDeformation_create_from_disp(const void* disp_ptr) {
	return cReg_NiftiImageData3DDeformation_create_from_disp(disp_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDeformation_get_inverse(const void* def_ptr, const void* floating_ptr) {
	return cReg_NiftiImageData3DDeformation_get_inverse(def_ptr, floating_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDisplacement_create_from_def(const void* def_ptr) {
	return cReg_NiftiImageData3DDisplacement_create_from_def(def_ptr);
}
EXPORTED_FUNCTION     void* mReg_Registration_process(void* ptr) {
	return cReg_Registration_process(ptr);
}
EXPORTED_FUNCTION     void* mReg_Registration_get_deformation_displacement_image(const void* ptr, const char *transform_type, const int idx) {
	return cReg_Registration_get_deformation_displacement_image(ptr, transform_type, idx);
}
EXPORTED_FUNCTION     void* mReg_Registration_add_floating(const void* ptr, const void *im_ptr) {
	return cReg_Registration_add_floating(ptr, im_ptr);
}
EXPORTED_FUNCTION     void* mReg_Registration_clear_floatings(const void* ptr) {
	return cReg_Registration_clear_floatings(ptr);
}
EXPORTED_FUNCTION     void* mReg_Registration_get_output(const void* ptr,const int idx) {
	return cReg_Registration_get_output(ptr, idx);
}
EXPORTED_FUNCTION     void* mReg_Registration_set_reference_image_filename(const void* ptr, const char* filename) {
	return cReg_Registration_set_reference_image_filename(ptr, filename);
}
EXPORTED_FUNCTION     void* mReg_Registration_set_floating_image_filename(const void* ptr, const char* filename) {
	return cReg_Registration_set_floating_image_filename(ptr, filename);
}
EXPORTED_FUNCTION     void* mReg_Registration_add_floating_image_filename(const void* ptr, const char* filename) {
	return cReg_Registration_add_floating_image_filename(ptr, filename);
}
EXPORTED_FUNCTION     void* mReg_NiftyRegistration_set_parameter(const void* ptr, const char* par, const char* arg1, const char* arg2) {
	return cReg_NiftyRegistration_set_parameter(ptr, par, arg1, arg2);
}
EXPORTED_FUNCTION     void* mReg_NiftyRegistration_print_all_wrapped_methods(const char* name) {
	return cReg_NiftyRegistration_print_all_wrapped_methods(name);
}
EXPORTED_FUNCTION     void* mReg_NiftyAladin_get_TM(const void* ptr, const char* dir) {
	return cReg_NiftyAladin_get_TM(ptr, dir);
}
EXPORTED_FUNCTION     void* mReg_SPMRegistration_get_TM(const void* ptr, const char* dir, const int idx) {
	return cReg_SPMRegistration_get_TM(ptr, dir, idx);
}
EXPORTED_FUNCTION     void* mReg_NiftyResample_add_transformation(void* self, const void* trans, const char* type) {
	return cReg_NiftyResample_add_transformation(self, trans, type);
}
EXPORTED_FUNCTION     void* mReg_NiftyResample_clear_transformations(void* self) {
	return cReg_NiftyResample_clear_transformations(self);
}
EXPORTED_FUNCTION     void* mReg_NiftyResample_process(void* ptr) {
	return cReg_NiftyResample_process(ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftyResample_forward(const void *output_ptr, const void * const input_ptr, const void *resampler_ptr) {
	return cReg_NiftyResample_forward(output_ptr, input_ptr, resampler_ptr);
}
EXPORTED_FUNCTION     void* mReg_NiftyResample_adjoint(const void *output_ptr, const void * const input_ptr, const void *resampler_ptr) {
	return cReg_NiftyResample_adjoint(output_ptr, input_ptr, resampler_ptr);
}
EXPORTED_FUNCTION     void* mReg_ImageWeightedMean_add_image(void* ptr, const void* obj, const float weight) {
	return cReg_ImageWeightedMean_add_image(ptr, obj, weight);
}
EXPORTED_FUNCTION     void* mReg_ImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight) {
	return cReg_ImageWeightedMean_add_image_filename(ptr, filename, weight);
}
EXPORTED_FUNCTION     void* mReg_ImageWeightedMean_process(void* ptr) {
	return cReg_ImageWeightedMean_process(ptr);
}
EXPORTED_FUNCTION     void* mReg_Transformation_get_as_deformation_field(const void* ptr, const char* name, const void* ref) {
	return cReg_Transformation_get_as_deformation_field(ptr, name, ref);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_construct_from_TM(PTR_FLOAT ptr_TM) {
	return cReg_AffineTransformation_construct_from_TM(ptr_TM);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_construct_from_trans_and_quaternion(PTR_FLOAT trans_ptr, const void* quat_ptr) {
	return cReg_AffineTransformation_construct_from_trans_and_quaternion(trans_ptr, quat_ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_construct_from_trans_and_euler(PTR_FLOAT trans_ptr, PTR_FLOAT euler_ptr) {
	return cReg_AffineTransformation_construct_from_trans_and_euler(trans_ptr, euler_ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_deep_copy(const void* ptr) {
	return cReg_AffineTransformation_deep_copy(ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_write(const void* ptr, const char* filename) {
	return cReg_AffineTransformation_write(ptr, filename);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_as_array(const void* ptr, PTR_FLOAT ptr_TM) {
	return cReg_AffineTransformation_as_array(ptr, ptr_TM);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_get_identity() {
	return cReg_AffineTransformation_get_identity();
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_get_inverse(const void* ptr) {
	return cReg_AffineTransformation_get_inverse(ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_get_Euler_angles(const void* ptr, PTR_FLOAT Euler) {
	return cReg_AffineTransformation_get_Euler_angles(ptr, Euler);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_get_quaternion(const void* ptr) {
	return cReg_AffineTransformation_get_quaternion(ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_mul(const void* mat1_ptr, const void* mat2_ptr) {
	return cReg_AffineTransformation_mul(mat1_ptr, mat2_ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_equal(const void* mat1_ptr, const void* mat2_ptr) {
	return cReg_AffineTransformation_equal(mat1_ptr, mat2_ptr);
}
EXPORTED_FUNCTION     void* mReg_AffineTransformation_get_average(const void *handle_vector_ptr) {
	return cReg_AffineTransformation_get_average(handle_vector_ptr);
}
EXPORTED_FUNCTION     void* mReg_Quaternion_construct_from_array(PTR_FLOAT arr) {
	return cReg_Quaternion_construct_from_array(arr);
}
EXPORTED_FUNCTION     void* mReg_Quaternion_construct_from_AffineTransformation(const void* ptr) {
	return cReg_Quaternion_construct_from_AffineTransformation(ptr);
}
EXPORTED_FUNCTION     void* mReg_Quaternion_get_average(const void *handle_vector_ptr) {
	return cReg_Quaternion_get_average(handle_vector_ptr);
}
EXPORTED_FUNCTION     void* mReg_Quaternion_as_array(const void* ptr, PTR_FLOAT arr) {
	return cReg_Quaternion_as_array(ptr, arr);
}
#ifndef CREG_FOR_MATLAB
}
#endif
void* newMexPrinter();
void* deleteMexPrinter(void* ptr);
EXPORTED_FUNCTION void* mNewMexPrinter() {
  return newMexPrinter();
}
EXPORTED_FUNCTION void* mDeleteMexPrinter(void* ptr) {
  return deleteMexPrinter(ptr);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
