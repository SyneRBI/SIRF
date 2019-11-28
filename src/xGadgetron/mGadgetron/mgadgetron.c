/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
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
#define CGADGETRON_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

#include <mex.h>
#include "matrix.h"
#include "sirf/Gadgetron/cgadgetron.h"

#ifndef CGADGETRON_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t
 extern "C" {
#else
#define PTR_INT int*
#define PTR_FLOAT float*
#define PTR_DOUBLE double*
#endif
EXPORTED_FUNCTION  void* mParameter(void* ptr, const char* obj, const char* name) {
	return parameter(ptr, obj, name);
}
EXPORTED_FUNCTION 	void* mSetParameter (void* ptr, const char* obj, const char* par, const void* val) {
	return setParameter (ptr, obj, par, val);
}
EXPORTED_FUNCTION 	void* mGT_newObject(const char* name) {
	return cGT_newObject(name);
}
EXPORTED_FUNCTION 	void* mGT_parameter(void* ptr, const char* obj, const char* name) {
	return cGT_parameter(ptr, obj, name);
}
EXPORTED_FUNCTION 	void* mGT_setParameter (void* ptr, const char* obj, const char* par, const void* val) {
	return cGT_setParameter (ptr, obj, par, val);
}
EXPORTED_FUNCTION 	void*	mGT_computeCoilImages(void* ptr_cis, void* ptr_acqs) {
	return cGT_computeCoilImages(ptr_cis, ptr_acqs);
}
EXPORTED_FUNCTION 	void*	mGT_computeCSMsFromCIs(void* ptr_csms, void* ptr_cis) {
	return cGT_computeCSMsFromCIs(ptr_csms, ptr_cis);
}
EXPORTED_FUNCTION 	void* mGT_CoilSensitivities(const char* file) {
	return cGT_CoilSensitivities(file);
}
EXPORTED_FUNCTION 	void* mGT_computeCoilSensitivities(void* ptr_csms, void* ptr_acqs) {
	return cGT_computeCoilSensitivities(ptr_csms, ptr_acqs);
}
EXPORTED_FUNCTION 	void* mGT_appendCSM (void* ptr_csms, int nx, int ny, int nz, int nc,  PTR_FLOAT ptr_re, PTR_FLOAT ptr_im) {
	return cGT_appendCSM (ptr_csms, nx, ny, nz, nc, ptr_re, ptr_im);
}
EXPORTED_FUNCTION 	void mGT_getCoilDataDimensions (void* ptr_csms, int csm_num, PTR_INT ptr_dim) {
	cGT_getCoilDataDimensions (ptr_csms, csm_num, ptr_dim);
}
EXPORTED_FUNCTION 	void mGT_getCoilData (void* ptr_csms, int csm_num, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im) {
	cGT_getCoilData (ptr_csms, csm_num, ptr_re, ptr_im);
}
EXPORTED_FUNCTION 	void* mGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs) {
	return cGT_AcquisitionModel(ptr_acqs, ptr_imgs);
}
EXPORTED_FUNCTION 	void* mGT_setUpAcquisitionModel (void* ptr_am, const void* ptr_acqs, const void* ptr_imgs) {
	return cGT_setUpAcquisitionModel (ptr_am, ptr_acqs, ptr_imgs);
}
EXPORTED_FUNCTION 	void* mGT_setAcquisitionModelParameter (void* ptr_am, const char* name, const void* ptr) {
	return cGT_setAcquisitionModelParameter (ptr_am, name, ptr);
}
EXPORTED_FUNCTION 	void* mGT_setCSMs(void* ptr_am, const void* ptr_csms) {
	return cGT_setCSMs(ptr_am, ptr_csms);
}
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs) {
	return cGT_AcquisitionModelForward(ptr_am, ptr_imgs);
}
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs) {
	return cGT_AcquisitionModelBackward(ptr_am, ptr_acqs);
}
EXPORTED_FUNCTION 	void* mGT_setAcquisitionDataStorageScheme(const char* scheme) {
	return cGT_setAcquisitionDataStorageScheme(scheme);
}
EXPORTED_FUNCTION 	void* mGT_getAcquisitionDataStorageScheme() {
	return cGT_getAcquisitionDataStorageScheme();
}
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFromFile(const char* file) {
	return cGT_ISMRMRDAcquisitionsFromFile(file);
}
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFile(const char* file) {
	return cGT_ISMRMRDAcquisitionsFile(file);
}
EXPORTED_FUNCTION 	void* mGT_processAcquisitions(void* ptr_proc, void* ptr_input) {
	return cGT_processAcquisitions(ptr_proc, ptr_input);
}
EXPORTED_FUNCTION 	void* mGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num) {
	return cGT_acquisitionFromContainer(ptr_acqs, acq_num);
}
EXPORTED_FUNCTION 	void* mGT_cloneAcquisitions(void* ptr_input) {
	return cGT_cloneAcquisitions(ptr_input);
}
EXPORTED_FUNCTION 	void* mGT_sortAcquisitions(void* ptr_acqs) {
	return cGT_sortAcquisitions(ptr_acqs);
}
EXPORTED_FUNCTION 	void* mGT_getAcquisitionDataDimensions(void* ptr_acqs, PTR_INT ptr_dim) {
	return cGT_getAcquisitionDataDimensions(ptr_acqs, ptr_dim);
}
EXPORTED_FUNCTION 	void* mGT_writeAcquisitions(void* ptr_acqs, const char* filename) {
	return cGT_writeAcquisitions(ptr_acqs, filename);
}
EXPORTED_FUNCTION 	void* mGT_fillAcquisitionData(void* ptr_acqs, PTR_FLOAT ptr_z, int all) {
	return cGT_fillAcquisitionData(ptr_acqs, ptr_z, all);
}
EXPORTED_FUNCTION 	void* mGT_fillAcquisitionDataFromAcquisitionData(void* ptr_dst, void* ptr_src) {
	return cGT_fillAcquisitionDataFromAcquisitionData(ptr_dst, ptr_src);
}
EXPORTED_FUNCTION 	void* mGT_acquisitionDataAsArray(void* ptr_acqs, PTR_FLOAT ptr_z, int all) {
	return cGT_acquisitionDataAsArray(ptr_acqs, ptr_z, all);
}
EXPORTED_FUNCTION 	void* mGT_reconstructImages(void* ptr_recon, void* ptr_input) {
	return cGT_reconstructImages(ptr_recon, ptr_input);
}
EXPORTED_FUNCTION 	void* mGT_reconstructedImages(void* ptr_recon) {
	return cGT_reconstructedImages(ptr_recon);
}
EXPORTED_FUNCTION 	void*	mGT_readImages(const char* file) {
	return cGT_readImages(file);
}
EXPORTED_FUNCTION 	void* mGT_processImages(void* ptr_proc, void* ptr_input) {
	return cGT_processImages(ptr_proc, ptr_input);
}
EXPORTED_FUNCTION 	void* mGT_selectImages (void* ptr_input, const char* attr, const char* target) {
	return cGT_selectImages (ptr_input, attr, target);
}
EXPORTED_FUNCTION 	void* mGT_writeImages (void* ptr_imgs, const char* out_file, const char* out_group) {
	return cGT_writeImages (ptr_imgs, out_file, out_group);
}
EXPORTED_FUNCTION 	void* mGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num) {
	return cGT_imageWrapFromContainer(ptr_imgs, img_num);
}
EXPORTED_FUNCTION 	void* mGT_imageDataType(const void* ptr_x, int im_num) {
	return cGT_imageDataType(ptr_x, im_num);
}
EXPORTED_FUNCTION 	void mGT_getImageDim(void* ptr_img, PTR_INT ptr_dim) {
	cGT_getImageDim(ptr_img, ptr_dim);
}
EXPORTED_FUNCTION 	void* mGT_imageType(const void* ptr_img) {
	return cGT_imageType(ptr_img);
}
EXPORTED_FUNCTION 	void* mGT_getImageDataAsFloatArray(void* ptr_imgs, PTR_FLOAT ptr_data) {
	return cGT_getImageDataAsFloatArray(ptr_imgs, ptr_data);
}
EXPORTED_FUNCTION 	void* mGT_setImageDataFromFloatArray(void* ptr_imgs, PTR_FLOAT ptr_data) {
	return cGT_setImageDataFromFloatArray(ptr_imgs, ptr_data);
}
EXPORTED_FUNCTION 	void* mGT_getImageDataAsCmplxArray(void* ptr_imgs, PTR_FLOAT ptr_z) {
	return cGT_getImageDataAsCmplxArray(ptr_imgs, ptr_z);
}
EXPORTED_FUNCTION 	void* mGT_setImageDataFromCmplxArray(void* ptr_imgs, PTR_FLOAT ptr_z) {
	return cGT_setImageDataFromCmplxArray(ptr_imgs, ptr_z);
}
EXPORTED_FUNCTION     void* mGT_print_header(const void* ptr_imgs, const int im_idx) {
	return cGT_print_header(ptr_imgs, im_idx);
}
EXPORTED_FUNCTION 	void* mGT_setHost(void* ptr_gc, const char* host) {
	return cGT_setHost(ptr_gc, host);
}
EXPORTED_FUNCTION 	void* mGT_setPort(void* ptr_gc, const char* port) {
	return cGT_setPort(ptr_gc, port);
}
EXPORTED_FUNCTION 	void* mGT_addReader(void* ptr_gc, const char* id, const void* ptr_r) {
	return cGT_addReader(ptr_gc, id, ptr_r);
}
EXPORTED_FUNCTION 	void* mGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r) {
	return cGT_addWriter(ptr_gc, id, ptr_r);
}
EXPORTED_FUNCTION 	void* mGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r) {
	return cGT_addGadget(ptr_gc, id, ptr_r);
}
EXPORTED_FUNCTION 	void* mGT_setGadgetProperty(void* ptr_g, const char* prop, const char* val) {
	return cGT_setGadgetProperty(ptr_g, prop, val);
}
EXPORTED_FUNCTION 	void* mGT_setGadgetProperties(void* ptr_g, const char* props) {
	return cGT_setGadgetProperties(ptr_g, props);
}
EXPORTED_FUNCTION 	void* mGT_configGadgetChain(void* ptr_con, void* ptr_gc) {
	return cGT_configGadgetChain(ptr_con, ptr_gc);
}
EXPORTED_FUNCTION 	void* mGT_registerImagesReceiver(void* ptr_con, void* ptr_img) {
	return cGT_registerImagesReceiver(ptr_con, ptr_img);
}
EXPORTED_FUNCTION 	void* mGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms) {
	return cGT_setConnectionTimeout(ptr_con, timeout_ms);
}
EXPORTED_FUNCTION 	void* mGT_connect(void* ptr_con, const char* host, const char* port) {
	return cGT_connect(ptr_con, host, port);
}
EXPORTED_FUNCTION 	void* mGT_sendConfigScript(void* ptr_con, const char* config) {
	return cGT_sendConfigScript(ptr_con, config);
}
EXPORTED_FUNCTION 	void* mGT_sendConfigFile(void* ptr_con, const char* file) {
	return cGT_sendConfigFile(ptr_con, file);
}
EXPORTED_FUNCTION 	void* mGT_sendParameters(void* ptr_con, const void* par) {
	return cGT_sendParameters(ptr_con, par);
}
EXPORTED_FUNCTION 	void* mGT_sendParametersString(void* ptr_con, const char* par) {
	return cGT_sendParametersString(ptr_con, par);
}
EXPORTED_FUNCTION 	void* mGT_sendAcquisitions(void* ptr_con, void* ptr_dat) {
	return cGT_sendAcquisitions(ptr_con, ptr_dat);
}
EXPORTED_FUNCTION 	void* mGT_sendImages(void* ptr_con, void* ptr_img) {
	return cGT_sendImages(ptr_con, ptr_img);
}
EXPORTED_FUNCTION 	void* mGT_disconnect(void* ptr_con) {
	return cGT_disconnect(ptr_con);
}
#ifndef CGADGETRON_FOR_MATLAB
}
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
