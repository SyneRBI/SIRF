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
#define CSTIR_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#include <mex.h>
#include "matrix.h"
#include "cstir.h"

#ifndef CSTIR_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t
 extern "C" {
#else
#define PTR_INT int*
#define PTR_FLOAT float*
#define PTR_DOUBLE double*
#endif
EXPORTED_FUNCTION  void* mSTIR_newObject(const char* name) {
	return cSTIR_newObject(name);
}
EXPORTED_FUNCTION 	void* mSTIR_objectFromFile(const char* name, const char* filename) {
	return cSTIR_objectFromFile(name, filename);
}
EXPORTED_FUNCTION 	void* mSTIR_setParameter (void* ptr, const char* obj, const char* name, const void* value) {
	return cSTIR_setParameter (ptr, obj, name, value);
}
EXPORTED_FUNCTION 	void* mSTIR_parameter(const void* ptr, const char* obj, const char* name) {
	return cSTIR_parameter(ptr, obj, name);
}
EXPORTED_FUNCTION 	void* mSTIR_setListmodeToSinogramsInterval (void* ptr_acq, PTR_FLOAT ptr_data) {
	return cSTIR_setListmodeToSinogramsInterval (ptr_acq, ptr_data);
}
EXPORTED_FUNCTION 	void* mSTIR_setListmodeToSinogramsFlag (void* ptr_lm2s, const char* flag, int v) {
	return cSTIR_setListmodeToSinogramsFlag (ptr_lm2s, flag, v);
}
EXPORTED_FUNCTION 	void* mSTIR_setupListmodeToSinogramsConverter(void* ptr) {
	return cSTIR_setupListmodeToSinogramsConverter(ptr);
}
EXPORTED_FUNCTION 	void* mSTIR_convertListmodeToSinograms(void* ptr) {
	return cSTIR_convertListmodeToSinograms(ptr);
}
EXPORTED_FUNCTION 	void* mSTIR_computeRandoms(void* ptr) {
	return cSTIR_computeRandoms(ptr);
}
EXPORTED_FUNCTION 	void* mSTIR_applyImageDataProcessor(const void* ptr_p, void* ptr_d) {
	return cSTIR_applyImageDataProcessor(ptr_p, ptr_d);
}
EXPORTED_FUNCTION 	void* mSTIR_createPETAcquisitionSensitivityModel (const void* ptr_src, const char* src) {
	return cSTIR_createPETAcquisitionSensitivityModel (ptr_src, src);
}
EXPORTED_FUNCTION 	void* mSTIR_createPETAttenuationModel (const void* ptr_img, const void* ptr_am) {
	return cSTIR_createPETAttenuationModel (ptr_img, ptr_am);
}
EXPORTED_FUNCTION 	void* mSTIR_chainPETAcquisitionSensitivityModels (const void* ptr_first, const void* ptr_second) {
	return cSTIR_chainPETAcquisitionSensitivityModels (ptr_first, ptr_second);
}
EXPORTED_FUNCTION 	void* mSTIR_setupAcquisitionSensitivityModel(void* ptr_sm, void* ptr_ad) {
	return cSTIR_setupAcquisitionSensitivityModel(ptr_sm, ptr_ad);
}
EXPORTED_FUNCTION 	void* mSTIR_applyAcquisitionSensitivityModel (void* ptr_sm, void* ptr_ad, const char* job) {
	return cSTIR_applyAcquisitionSensitivityModel (ptr_sm, ptr_ad, job);
}
EXPORTED_FUNCTION 	void* mSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im) {
	return cSTIR_setupAcquisitionModel(ptr_am, ptr_dt, ptr_im);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelFwd(void* ptr_am, void* ptr_im,  int subset_num, int num_subsets) {
	return cSTIR_acquisitionModelFwd(ptr_am, ptr_im, subset_num, num_subsets);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad,  int subset_num, int num_subsets) {
	return cSTIR_acquisitionModelBwd(ptr_am, ptr_ad, subset_num, num_subsets);
}
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionsStorageScheme() {
	return cSTIR_getAcquisitionsStorageScheme();
}
EXPORTED_FUNCTION 	void* mSTIR_setAcquisitionsStorageScheme(const char* scheme) {
	return cSTIR_setAcquisitionsStorageScheme(scheme);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionsDataFromTemplate(void* ptr_t) {
	return cSTIR_acquisitionsDataFromTemplate(ptr_t);
}
EXPORTED_FUNCTION 	void* mSTIR_rebinnedAcquisitionData(void* ptr_t, const int num_segments_to_combine, const int num_views_to_combine, const int num_tang_poss_to_trim, const bool do_normalisation, const int max_in_segment_num_to_process ) {
	return cSTIR_rebinnedAcquisitionData(ptr_t, num_segments_to_combine, num_views_to_combine, num_tang_poss_to_trim, do_normalisation, max_in_segment_num_to_process);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionsDataFromScannerInfo (const char* scanner, int span, int max_ring_diff, int view_mash_factor) {
	return cSTIR_acquisitionsDataFromScannerInfo (scanner, span, max_ring_diff, view_mash_factor);
}
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionsDimensions(const void* ptr_acq, PTR_INT ptr_dim) {
	return cSTIR_getAcquisitionsDimensions(ptr_acq, ptr_dim);
}
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionsData(const void* ptr_acq, PTR_FLOAT ptr_data) {
	return cSTIR_getAcquisitionsData(ptr_acq, ptr_data);
}
EXPORTED_FUNCTION 	void* mSTIR_setAcquisitionsData(void* ptr_acq, PTR_FLOAT ptr_data) {
	return cSTIR_setAcquisitionsData(ptr_acq, ptr_data);
}
EXPORTED_FUNCTION 	void* mSTIR_fillAcquisitionsData(void* ptr_acq, float v) {
	return cSTIR_fillAcquisitionsData(ptr_acq, v);
}
EXPORTED_FUNCTION 	void* mSTIR_fillAcquisitionsDataFromAcquisitionsData (void* ptr_acq, const void * ptr_from) {
	return cSTIR_fillAcquisitionsDataFromAcquisitionsData (ptr_acq, ptr_from);
}
EXPORTED_FUNCTION 	void* mSTIR_writeAcquisitionData(void* ptr_acq, const char* filename) {
	return cSTIR_writeAcquisitionData(ptr_acq, filename);
}
EXPORTED_FUNCTION 	void* mSTIR_setupFBP2DReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_setupFBP2DReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_runFBP2DReconstruction(void* ptr_r) {
	return cSTIR_runFBP2DReconstruction(ptr_r);
}
EXPORTED_FUNCTION 	void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_setupReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_runReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_updateReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_setupObjectiveFunction(void* ptr_r, void* ptr_i) {
	return cSTIR_setupObjectiveFunction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void*	mSTIR_subsetSensitivity(void* ptr_f, int subset) {
	return cSTIR_subsetSensitivity(ptr_f, subset);
}
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionValue(void* ptr_f, void* ptr_i) {
	return cSTIR_objectiveFunctionValue(ptr_f, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionGradient (void* ptr_f, void* ptr_i, int subset) {
	return cSTIR_objectiveFunctionGradient (ptr_f, ptr_i, subset);
}
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionGradientNotDivided (void* ptr_f, void* ptr_i, int subset) {
	return cSTIR_objectiveFunctionGradientNotDivided (ptr_f, ptr_i, subset);
}
EXPORTED_FUNCTION 	void* mSTIR_setupPrior(void* ptr_p, void* ptr_i) {
	return cSTIR_setupPrior(ptr_p, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_priorGradient(void* ptr_p, void* ptr_i) {
	return cSTIR_priorGradient(ptr_p, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_getImageDimensions(const void* ptr, PTR_INT ptr_data) {
	return cSTIR_getImageDimensions(ptr, ptr_data);
}
EXPORTED_FUNCTION 	void* mSTIR_getImageVoxelSizes(const void* ptr_im, PTR_FLOAT ptr_vs) {
	return cSTIR_getImageVoxelSizes(ptr_im, ptr_vs);
}
EXPORTED_FUNCTION 	void* mSTIR_getImageData(const void* ptr, PTR_FLOAT ptr_data) {
	return cSTIR_getImageData(ptr, ptr_data);
}
EXPORTED_FUNCTION 	void* mSTIR_setImageData(const void* ptr_im, PTR_FLOAT ptr_data) {
	return cSTIR_setImageData(ptr_im, ptr_data);
}
EXPORTED_FUNCTION 	void* mSTIR_voxels3DF(int nx, int ny, int nz, float sx, float sy, float sz, float x, float y, float z) {
	return cSTIR_voxels3DF(nx, ny, nz, sx, sy, sz, x, y, z);
}
EXPORTED_FUNCTION 	void* mSTIR_imageFromVoxels(void* ptr_v) {
	return cSTIR_imageFromVoxels(ptr_v);
}
EXPORTED_FUNCTION 	void* mSTIR_imageFromImage(void* ptr_v) {
	return cSTIR_imageFromImage(ptr_v);
}
EXPORTED_FUNCTION 	void* mSTIR_imageFromAcquisitionData(void* ptr_ad) {
	return cSTIR_imageFromAcquisitionData(ptr_ad);
}
EXPORTED_FUNCTION 	void* mSTIR_imageFromAcquisitionDataAndNxNy(void* ptr_ad, int nx, int ny) {
	return cSTIR_imageFromAcquisitionDataAndNxNy(ptr_ad, nx, ny);
}
EXPORTED_FUNCTION 	void* mSTIR_fillImage(void* ptr_i, float v) {
	return cSTIR_fillImage(ptr_i, v);
}
EXPORTED_FUNCTION 	void* mSTIR_addShape(void* ptr_i, void* ptr_s, float v) {
	return cSTIR_addShape(ptr_i, ptr_s, v);
}
EXPORTED_FUNCTION 	void* mSTIR_writeImage(void* ptr_i, const char* filename) {
	return cSTIR_writeImage(ptr_i, filename);
}
EXPORTED_FUNCTION 	void* mSTIR_norm(const void* ptr_x) {
	return cSTIR_norm(ptr_x);
}
EXPORTED_FUNCTION 	void*	mSTIR_dot(const void* ptr_x, const void* ptr_y) {
	return cSTIR_dot(ptr_x, ptr_y);
}
EXPORTED_FUNCTION 	void* mSTIR_mult(float a, const void* ptr_x) {
	return cSTIR_mult(a, ptr_x);
}
EXPORTED_FUNCTION 	void* mSTIR_axpby(float a, const void* ptr_x, float b, const void* ptr_y) {
	return cSTIR_axpby(a, ptr_x, b, ptr_y);
}
EXPORTED_FUNCTION 	void* mNewTextPrinter(const char* stream) {
	return newTextPrinter(stream);
}
EXPORTED_FUNCTION 	void* mNewTextWriter(const char* stream) {
	return newTextWriter(stream);
}
EXPORTED_FUNCTION 	void mOpenChannel(int channel, void* ptr_w) {
	openChannel(channel, ptr_w);
}
EXPORTED_FUNCTION 	void mCloseChannel(int channel, void* ptr_w) {
	closeChannel(channel, ptr_w);
}
EXPORTED_FUNCTION 	void* mDeleteTextPrinter(void* ptr) {
	return deleteTextPrinter(ptr);
}
EXPORTED_FUNCTION 	void* mDeleteTextWriter(void* ptr_w) {
	return deleteTextWriter(ptr_w);
}
#ifndef CSTIR_FOR_MATLAB
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
