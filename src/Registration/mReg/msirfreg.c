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
#define CSIRFREG_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#include <mex.h>
#include "matrix.h"
#include "csirfreg.h"

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
EXPORTED_FUNCTION  void* mSIRFReg_newObject(const char* name) {
	return cSIRFReg_newObject(name);
}
EXPORTED_FUNCTION 	void* mSIRFReg_objectFromFile(const char* name, const char* filename) {
	return cSIRFReg_objectFromFile(name, filename);
}
EXPORTED_FUNCTION 	void* mSIRFReg_setParameter (void* ptr, const char* obj, const char* name, const void* value) {
	return cSIRFReg_setParameter (ptr, obj, name, value);
}
EXPORTED_FUNCTION 	void* mSIRFReg_parameter(const void* ptr, const char* obj, const char* name) {
	return cSIRFReg_parameter(ptr, obj, name);
}
EXPORTED_FUNCTION     void* mSIRFReg_do_nifti_images_match(const void* im1, const void* im2, const float accuracy_percentage_of_max) {
	return cSIRFReg_do_nifti_images_match(im1, im2, accuracy_percentage_of_max);
}
EXPORTED_FUNCTION     void* mSIRFReg_dump_nifti_info_filename(const char* filename) {
	return cSIRFReg_dump_nifti_info_filename(filename);
}
EXPORTED_FUNCTION     void* mSIRFReg_dump_nifti_info_im1(const void* im1) {
	return cSIRFReg_dump_nifti_info_im1(im1);
}
EXPORTED_FUNCTION     void* mSIRFReg_dump_nifti_info_im2(const void* im1, const void* im2) {
	return cSIRFReg_dump_nifti_info_im2(im1, im2);
}
EXPORTED_FUNCTION     void* mSIRFReg_dump_nifti_info_im3(const void* im1, const void* im2, const void* im3) {
	return cSIRFReg_dump_nifti_info_im3(im1, im2, im3);
}
EXPORTED_FUNCTION     void* mSIRFReg_dump_nifti_info_im4(const void* im1, const void* im2, const void* im3, const void* im4) {
	return cSIRFReg_dump_nifti_info_im4(im1, im2, im3, im4);
}
EXPORTED_FUNCTION     void* mSIRFReg_dump_nifti_info_im5(const void* im1, const void* im2, const void* im3, const void* im4, const void* im5) {
	return cSIRFReg_dump_nifti_info_im5(im1, im2, im3, im4, im5);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFReg_open_TM(const char* filename, PTR_FLOAT ptr_TM) {
	return cSIRFReg_SIRFReg_open_TM(filename, ptr_TM);
}
EXPORTED_FUNCTION     void* mSIRFReg_compose_transformations_into_single_deformation2(const void* im, const void* trans1, const void* trans2) {
	return cSIRFReg_compose_transformations_into_single_deformation2(im, trans1, trans2);
}
EXPORTED_FUNCTION     void* mSIRFReg_compose_transformations_into_single_deformation3(const void* im, const void* trans1, const void* trans2, const void* trans3) {
	return cSIRFReg_compose_transformations_into_single_deformation3(im, trans1, trans2, trans3);
}
EXPORTED_FUNCTION     void* mSIRFReg_compose_transformations_into_single_deformation4(const void* im, const void* trans1, const void* trans2, const void* trans3, const void* trans4) {
	return cSIRFReg_compose_transformations_into_single_deformation4(im, trans1, trans2, trans3, trans4);
}
EXPORTED_FUNCTION     void* mSIRFReg_compose_transformations_into_single_deformation5(const void* im, const void* trans1, const void* trans2, const void* trans3, const void* trans4, const void* trans5) {
	return cSIRFReg_compose_transformations_into_single_deformation5(im, trans1, trans2, trans3, trans4, trans5);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_from_PETImageData(void* ptr) {
	return cSIRFReg_SIRFImageData_from_PETImageData(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_save_to_file(const void* ptr, const char* filename) {
	return cSIRFReg_SIRFImageData_save_to_file(ptr, filename);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_copy_data_to(const void* ptr, const void* obj) {
	return cSIRFReg_SIRFImageData_copy_data_to(ptr, obj);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_fill(const void* ptr, const float val) {
	return cSIRFReg_SIRFImageData_fill(ptr, val);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_deep_copy(const void* ptr) {
	return cSIRFReg_SIRFImageData_deep_copy(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_get_dimensions(const void* ptr, PTR_INT ptr_dim) {
	return cSIRFReg_SIRFImageData_get_dimensions(ptr, ptr_dim);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_get_data(const void* ptr, PTR_FLOAT ptr_data) {
	return cSIRFReg_SIRFImageData_get_data(ptr, ptr_data);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageData_maths(const void* ptr, const void* obj, const int maths_type) {
	return cSIRFReg_SIRFImageData_maths(ptr, obj, maths_type);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageDataDeformation_save_to_file(const void* ptr, const char* filename, const bool split_xyz) {
	return cSIRFReg_SIRFImageDataDeformation_save_to_file(ptr, filename, split_xyz);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageDataDeformation_create_from_3D_image(void* ptr, const void* obj) {
	return cSIRFReg_SIRFImageDataDeformation_create_from_3D_image(ptr, obj);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFImageDataDeformation_deep_copy(const void* ptr) {
	return cSIRFReg_SIRFImageDataDeformation_deep_copy(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFReg_update(void* ptr) {
	return cSIRFReg_SIRFReg_update(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFReg_get_deformation_displacement_image(const void* ptr, const char* type) {
	return cSIRFReg_SIRFReg_get_deformation_displacement_image(ptr, type);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(const void* ptr, const char* filename, const char* dir) {
	return cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(ptr, filename, dir);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFReg_get_TM(const void* ptr, PTR_FLOAT ptr_TM, const char* dir) {
	return cSIRFReg_SIRFReg_get_TM(ptr, ptr_TM, dir);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegNiftyResample_add_transformation(void* self, const void* trans, const char *type) {
	return cSIRFReg_SIRFRegNiftyResample_add_transformation(self, trans, type);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegNiftyResample_update(void* ptr) {
	return cSIRFReg_SIRFRegNiftyResample_update(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegImageWeightedMean3D_add_image(void* ptr, const void* obj, const float weight) {
	return cSIRFReg_SIRFRegImageWeightedMean3D_add_image(ptr, obj, weight);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegImageWeightedMean3D_add_image_filename(void* ptr, const char* filename, const float weight) {
	return cSIRFReg_SIRFRegImageWeightedMean3D_add_image_filename(ptr, filename, weight);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegImageWeightedMean3D_update(void* ptr) {
	return cSIRFReg_SIRFRegImageWeightedMean3D_update(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegImageWeightedMean4D_add_image(void* ptr, const void* obj, const float weight) {
	return cSIRFReg_SIRFRegImageWeightedMean4D_add_image(ptr, obj, weight);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegImageWeightedMean4D_add_image_filename(void* ptr, const char* filename, const float weight) {
	return cSIRFReg_SIRFRegImageWeightedMean4D_add_image_filename(ptr, filename, weight);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegImageWeightedMean4D_update(void* ptr) {
	return cSIRFReg_SIRFRegImageWeightedMean4D_update(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegTransformation_get_as_deformation_field(const void* ptr, const void* ref) {
	return cSIRFReg_SIRFRegTransformation_get_as_deformation_field(ptr, ref);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegTransformationAffine_construct_from_TM(PTR_FLOAT ptr_TM) {
	return cSIRFReg_SIRFRegTransformationAffine_construct_from_TM(ptr_TM);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegTransformationDisplacement_construct_from_SIRFImageDataDeformation(const void* ptr) {
	return cSIRFReg_SIRFRegTransformationDisplacement_construct_from_SIRFImageDataDeformation(ptr);
}
EXPORTED_FUNCTION     void* mSIRFReg_SIRFRegTransformationDeformation_construct_from_SIRFImageDataDeformation(const void* ptr) {
	return cSIRFReg_SIRFRegTransformationDeformation_construct_from_SIRFImageDataDeformation(ptr);
}
#ifndef CSIRFREG_FOR_MATLAB
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
