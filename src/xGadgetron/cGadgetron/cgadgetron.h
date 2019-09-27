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

#ifndef GADGETRON_C_INTERFACE
#define GADGETRON_C_INTERFACE

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

	// Unified parameter exchange methods
	void* parameter(void* ptr, const char* obj, const char* name);
	void* setParameter
		(void* ptr, const char* obj, const char* par, const void* val);
	// common Object methods
	void* cGT_newObject(const char* name);
	void* cGT_parameter(void* ptr, const char* obj, const char* name);
	void* cGT_setParameter
		(void* ptr, const char* obj, const char* par, const void* val);

	// coil data methods
	void*	cGT_computeCoilImages(void* ptr_cis, void* ptr_acqs);
	void*	cGT_computeCSMsFromCIs(void* ptr_csms, void* ptr_cis);
	void* cGT_CoilSensitivities(const char* file);
	void* cGT_computeCoilSensitivities(void* ptr_csms, void* ptr_acqs);
	void* cGT_appendCSM
		(void* ptr_csms, int nx, int ny, int nz, int nc, 
		PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);
	void cGT_getCoilDataDimensions
		(void* ptr_csms, int csm_num, PTR_INT ptr_dim);
	void cGT_getCoilData
		(void* ptr_csms, int csm_num, PTR_FLOAT ptr_re, PTR_FLOAT ptr_im);

	// acquisition model methods
	void* cGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs);
	void* cGT_setUpAcquisitionModel
		(void* ptr_am, const void* ptr_acqs, const void* ptr_imgs);
	void* cGT_setAcquisitionModelParameter
		(void* ptr_am, const char* name, const void* ptr);
	void* cGT_setCSMs(void* ptr_am, const void* ptr_csms);
	void* cGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
	void* cGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs);

	// acquisition data methods
	void* cGT_setAcquisitionDataStorageScheme(const char* scheme);
	void* cGT_getAcquisitionDataStorageScheme();
	void* cGT_ISMRMRDAcquisitionsFromFile(const char* file);
	void* cGT_ISMRMRDAcquisitionsFile(const char* file);
	void* cGT_processAcquisitions(void* ptr_proc, void* ptr_input);
	void* cGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num);
	void* cGT_cloneAcquisitions(void* ptr_input);
	void* cGT_sortAcquisitions(void* ptr_acqs);
	void* cGT_getAcquisitionDataDimensions(void* ptr_acqs, PTR_INT ptr_dim);
	void* cGT_writeAcquisitions(void* ptr_acqs, const char* filename);
	void* cGT_fillAcquisitionData(void* ptr_acqs, PTR_FLOAT ptr_z, int all);
	void* cGT_fillAcquisitionDataFromAcquisitionData(void* ptr_dst, void* ptr_src);
	void* cGT_acquisitionDataAsArray(void* ptr_acqs, PTR_FLOAT ptr_z, int all);

	// image methods
	void* cGT_reconstructImages(void* ptr_recon, void* ptr_input);
	void* cGT_reconstructedImages(void* ptr_recon);
	void*	cGT_readImages(const char* file);
	void* cGT_processImages(void* ptr_proc, void* ptr_input);
	void* cGT_selectImages
		(void* ptr_input, const char* attr, const char* target);
	void* cGT_writeImages
		(void* ptr_imgs, const char* out_file, const char* out_group);
	void* cGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num);
	void* cGT_imageDataType(const void* ptr_x, int im_num);
	void cGT_getImageDim(void* ptr_img, PTR_INT ptr_dim);
	void* cGT_imageType(const void* ptr_img);
	void* cGT_getImageDataAsFloatArray(void* ptr_imgs, PTR_FLOAT ptr_data);
	void* cGT_setImageDataFromFloatArray(void* ptr_imgs, PTR_FLOAT ptr_data);
	void* cGT_getImageDataAsCmplxArray(void* ptr_imgs, PTR_FLOAT ptr_z);
	void* cGT_setImageDataFromCmplxArray(void* ptr_imgs, PTR_FLOAT ptr_z);
    void* cGT_print_header(const void* ptr_imgs, const int im_idx);

	// gadget chain methods
	void* cGT_setHost(void* ptr_gc, const char* host);
	void* cGT_setPort(void* ptr_gc, const char* port);
	void* cGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_setGadgetProperty(void* ptr_g, const char* prop, const char* val);
	void* cGT_setGadgetProperties(void* ptr_g, const char* props);
	void* cGT_configGadgetChain(void* ptr_con, void* ptr_gc);
	void* cGT_registerImagesReceiver(void* ptr_con, void* ptr_img);

	// gadgetron client methods
	void* cGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms);
	void* cGT_connect(void* ptr_con, const char* host, const char* port);
	void* cGT_sendConfigScript(void* ptr_con, const char* config);
	void* cGT_sendConfigFile(void* ptr_con, const char* file);
	void* cGT_sendParameters(void* ptr_con, const void* par);
	void* cGT_sendParametersString(void* ptr_con, const char* par);
	void* cGT_sendAcquisitions(void* ptr_con, void* ptr_dat);
	void* cGT_sendImages(void* ptr_con, void* ptr_img);
	void* cGT_disconnect(void* ptr_con);

#ifndef CGADGETRON_FOR_MATLAB
}
#endif

#endif
