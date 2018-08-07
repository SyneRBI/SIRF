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
#ifndef CGADGETRON_TO_MATLAB_INTERFACE
#define CGADGETRON_TO_MATLAB_INTERFACE

#define CGADGETRON_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

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
EXPORTED_FUNCTION  void* mGT_newObject(const char* name);
EXPORTED_FUNCTION 	void* mGT_parameter(void* ptr, const char* obj, const char* name);
EXPORTED_FUNCTION 	void* mGT_setParameter (void* ptr, const char* obj, const char* par, const void* val);
EXPORTED_FUNCTION 	void*	mGT_computeCoilImages(void* ptr_cis, void* ptr_acqs);
EXPORTED_FUNCTION 	void*	mGT_computeCSMsFromCIs(void* ptr_csms, void* ptr_cis);
EXPORTED_FUNCTION 	void* mGT_CoilSensitivities(const char* file);
EXPORTED_FUNCTION 	void* mGT_computeCoilSensitivities(void* ptr_csms, void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_appendCSM (void* ptr_csms, int nx, int ny, int nz, int nc,  PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_setCSMs(void* ptr_am, const void* ptr_csms);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_setAcquisitionsStorageScheme(const char* scheme);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsStorageScheme();
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFromFile(const char* file);
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFile(const char* file);
EXPORTED_FUNCTION 	void* mGT_processAcquisitions(void* ptr_proc, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num);
EXPORTED_FUNCTION 	void* mGT_orderAcquisitions(void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsDimensions(void* ptr_acqs, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsFlags (void* ptr_acqs, unsigned int n, PTR_INT ptr_f);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsData (void* ptr_acqs, unsigned int slice, PTR_FLOAT ptr_r, PTR_FLOAT ptr_i);
EXPORTED_FUNCTION 	void* mGT_setAcquisitionsData (void* ptr_acqs, unsigned int na, unsigned int nc, unsigned int ns, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
EXPORTED_FUNCTION 	void*	mGT_writeAcquisitions(void* ptr_acqs, const char* filename);
EXPORTED_FUNCTION 	void* mGT_reconstructImages(void* ptr_recon, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_reconstructedImages(void* ptr_recon);
EXPORTED_FUNCTION 	void*	mGT_readImages(const char* file);
EXPORTED_FUNCTION 	void* mGT_processImages(void* ptr_proc, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_selectImages (void* ptr_input, const char* attr, const char* target);
EXPORTED_FUNCTION 	void* mGT_writeImages (void* ptr_imgs, const char* out_file, const char* out_group);
EXPORTED_FUNCTION 	void* mGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num);
EXPORTED_FUNCTION 	void* mGT_imageDataType(const void* ptr_x, int im_num);
EXPORTED_FUNCTION 	void mGT_getCoilDataDimensions (void* ptr_csms, int csm_num, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void mGT_getCoilData (void* ptr_csms, int csm_num, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
EXPORTED_FUNCTION 	void mGT_getCoilDataAbs(void* ptr_csms, int csm_num, PTR_FLOAT ptr);
EXPORTED_FUNCTION 	void mGT_getImageDim(void* ptr_img, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void* mGT_imageType(const void* ptr_img);
EXPORTED_FUNCTION 	void mGT_getImageDataAsFloatArray(void* ptr_img, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void mGT_getImageDataAsComplexArray (void* ptr_imgs, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
EXPORTED_FUNCTION 	void mGT_getImageDimensions(void* ptr_imgs, int img_num, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void mGT_getImagesDataAsFloatArray(void* ptr_imgs, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void mGT_getImagesDataAsComplexArray (void* ptr_imgs, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
EXPORTED_FUNCTION 	void* mGT_setComplexImagesData (void* ptr_imgs, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
EXPORTED_FUNCTION 	void* mGT_dataItems(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_norm(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_dot(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_axpby( float ar, float ai, const void* ptr_x, float br, float bi, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_multiply(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_divide(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_setGadgetProperty(void* ptr_g, const char* prop, const char* val);
EXPORTED_FUNCTION 	void* mGT_setGadgetProperties(void* ptr_g, const char* props);
EXPORTED_FUNCTION 	void* mGT_configGadgetChain(void* ptr_con, void* ptr_gc);
EXPORTED_FUNCTION 	void* mGT_registerImagesReceiver(void* ptr_con, void* ptr_img);
EXPORTED_FUNCTION 	void* mGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms);
EXPORTED_FUNCTION 	void* mGT_connect(void* ptr_con, const char* host, const char* port);
EXPORTED_FUNCTION 	void* mGT_sendConfigScript(void* ptr_con, const char* config);
EXPORTED_FUNCTION 	void* mGT_sendConfigFile(void* ptr_con, const char* file);
EXPORTED_FUNCTION 	void* mGT_sendParameters(void* ptr_con, const void* par);
EXPORTED_FUNCTION 	void* mGT_sendParametersString(void* ptr_con, const char* par);
EXPORTED_FUNCTION 	void* mGT_sendAcquisitions(void* ptr_con, void* ptr_dat);
EXPORTED_FUNCTION 	void* mGT_sendImages(void* ptr_con, void* ptr_img);
EXPORTED_FUNCTION 	void* mGT_disconnect(void* ptr_con);
#ifndef CGADGETRON_FOR_MATLAB
}
#endif

#endif
