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
#ifndef CSIRFREG_TO_MATLAB_INTERFACE
#define CSIRFREG_TO_MATLAB_INTERFACE

#define CSIRFREG_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

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
EXPORTED_FUNCTION  void* mReg_newObject(const char* name);
EXPORTED_FUNCTION     void* mReg_objectFromFile(const char* name, const char* filename);
EXPORTED_FUNCTION     void* mReg_setParameter(void* ptr, const char* obj, const char* name, const void* value);
EXPORTED_FUNCTION     void* mReg_parameter(const void* ptr, const char* obj, const char* name);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_print_headers(const int num_ims, const void* im1, const void* im2, const void* im3, const void* im4, const void* im5);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_write(const void* ptr, const char* filename, const int datatype);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_fill(const void* ptr, const float val);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_fill_arr(const void* ptr, PTR_FLOAT val);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_deep_copy(const void* copy_ptr, const void *orig_ptr);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_dimensions(const void* ptr, PTR_INT ptr_dim);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_data(const void* ptr, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_maths_num(const void *res_ptr, const void* im1_ptr, const float val, const int maths_type);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_equal(const void* im1_ptr, const void* im2_ptr);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_norm(const void* im1_ptr, const void* im2_ptr);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_get_original_datatype(const void* im_ptr);
EXPORTED_FUNCTION     void* mReg_NiftiImageData_crop(const void* im_ptr, PTR_INT min_index_ptr, PTR_INT max_index_ptr);
/* TODO UNCOMMENT WHEN GEOMETRICAL INFO IS IMPLEMENTED
    void* cReg_NiftiImageData3D_from_PETImageData(void* ptr);
    void* cReg_NiftiImageData3D_copy_data_to(const void* ptr, const void* obj);
*/
EXPORTED_FUNCTION  void* mReg_NiftiImageData3DTensor_write_split_xyz_components(const void* ptr, const char* filename, const int datatype);
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_create_from_3D_image(const void *ptr, const void* obj);
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr);
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DTensor_flip_component(const void *ptr, const int dim);
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDeformation_compose_single_deformation(const void* im, const int num_elements, const char* types, const void* trans1, const void* trans2, const void* trans3, const void* trans4, const void* trans5);
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDeformation_create_from_disp(const void* disp_ptr);
EXPORTED_FUNCTION     void* mReg_NiftiImageData3DDisplacement_create_from_def(const void* def_ptr);
EXPORTED_FUNCTION     void* mReg_SIRFReg_process(void* ptr);
EXPORTED_FUNCTION     void* mReg_SIRFReg_get_deformation_displacement_image(const void* ptr, const char *transform_type);
EXPORTED_FUNCTION     void* mReg_SIRFReg_set_parameter(const void* ptr, const char* par, const char* arg1, const char* arg2);
EXPORTED_FUNCTION     void* mReg_SIRFReg_get_TM(const void* ptr, const char* dir);
EXPORTED_FUNCTION     void* mReg_SIRFRegNiftyResample_add_transformation(void* self, const void* trans, const char* type);
EXPORTED_FUNCTION     void* mReg_SIRFRegNiftyResample_process(void* ptr);
EXPORTED_FUNCTION     void* mReg_SIRFRegImageWeightedMean_add_image(void* ptr, const void* obj, const float weight);
EXPORTED_FUNCTION     void* mReg_SIRFRegImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight);
EXPORTED_FUNCTION     void* mReg_SIRFRegImageWeightedMean_process(void* ptr);
EXPORTED_FUNCTION     void* mReg_SIRFRegTransformation_get_as_deformation_field(const void* ptr, const char* name, const void* ref);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_construct_from_TM(PTR_FLOAT ptr_TM);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_deep_copy(const void* ptr);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_write(const void* ptr, const char* filename);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_as_array(const void* ptr, PTR_FLOAT ptr_TM);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_get_identity();
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_get_inverse(const void* ptr);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_mul(const void* mat1_ptr, const void* mat2_ptr);
EXPORTED_FUNCTION     void* mReg_SIRFRegAffineTransformation_equal(const void* mat1_ptr, const void* mat2_ptr);
#ifndef CSIRFREG_FOR_MATLAB
}
#endif
EXPORTED_FUNCTION void* mNewMexPrinter();
EXPORTED_FUNCTION void* mDeleteMexPrinter(void* ptr);

#endif
