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
	void* cGT_newObject(const char* name);
	void* cGT_parameter(void* ptr, const char* obj, const char* name);
	void* cGT_setParameter
		(void* ptr, const char* obj, const char* par, const void* val);

	void*	cGT_computeCoilImages(void* ptr_cis, void* ptr_acqs);
	void*	cGT_computeCSMsFromCIs(void* ptr_csms, void* ptr_cis);
	void* cGT_CoilSensitivities(const char* file);
	void* cGT_computeCoilSensitivities(void* ptr_csms, void* ptr_acqs);
	void* cGT_appendCSM
		(void* ptr_csms, int nx, int ny, int nz, int nc, 
		PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);

	void* cGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs);
	void* cGT_setCSMs(void* ptr_am, const void* ptr_csms);
	void* cGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
	void* cGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs);

	void* cGT_ISMRMRDAcquisitionsFromFile(const char* file);
	void* cGT_ISMRMRDAcquisitionsFile(const char* file);
	void* cGT_processAcquisitions(void* ptr_proc, void* ptr_input);
	void* cGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num);
	void* cGT_orderAcquisitions(void* ptr_acqs);
	void* cGT_getAcquisitionsDimensions(void* ptr_acqs, PTR_INT ptr_dim);
	void* cGT_getAcquisitionsFlags(void* ptr_acqs, unsigned int n, PTR_INT ptr_f);
	void* cGT_getAcquisitionsData
		(void* ptr_acqs, unsigned int slice, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
	void* cGT_setAcquisitionsData
		(void* ptr_acqs, unsigned int na, unsigned int nc, unsigned int ns,
		PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);

	void* cGT_reconstructImages(void* ptr_recon, void* ptr_input);
	void* cGT_reconstructedImages(void* ptr_recon);
	void* cGT_processImages(void* ptr_proc, void* ptr_input);
	void* cGT_selectImages(void* ptr_input, const char* attr, const char* target);
	void cGT_setImageToRealConversion(void* ptr_imgs, int type);
	void* cGT_imagesCopy(const void* ptr_imgs);
	void* cGT_writeImages
		(void* ptr_imgs, const char* out_file, const char* out_group);
	void* cGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num);
	void* cGT_imageTypes(const void* ptr_x);
	void* cGT_imageDataType(const void* ptr_x, int im_num);

	void cGT_getCoilDataDimensions(void* ptr_csms, int csm_num, PTR_INT ptr_dim);
	void cGT_getCoilData
		(void* ptr_csms, int csm_num, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
	void cGT_getCoilDataAbs(void* ptr_csms, int csm_num, PTR_DOUBLE ptr);
	void cGT_getImageDimensions(void* ptr_imgs, int img_num, PTR_INT ptr_dim);
	void cGT_getImageDataAsDoubleArray
		(void* ptr_imgs, int img_num, PTR_DOUBLE ptr_data);
	void cGT_getImageDataAsComplexArray
		(void* ptr_imgs, int img_num, PTR_DOUBLE ptr_data);
	void cGT_getImagesDataAsDoubleArray(void* ptr_imgs, PTR_DOUBLE ptr_data);
	void cGT_getImagesDataAsComplexArray
		(void* ptr_imgs, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
	void* cGT_setComplexImagesData(void* ptr_imgs, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);

	void* cGT_dataItems(const void* ptr_x);
	void* cGT_norm(const void* ptr_x);
	void* cGT_dot(const void* ptr_x, const void* ptr_y);
	void* cGT_axpby(
		double ar, double ai, const void* ptr_x,
		double br, double bi, const void* ptr_y);

	void* cGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_setGadgetProperty(void* ptr_g, const char* prop, const char* value);
	void* cGT_setGadgetProperties(void* ptr_g, const char* props);
	void* cGT_configGadgetChain(void* ptr_con, void* ptr_gc);
	void* cGT_registerImagesReceiver(void* ptr_con, void* ptr_img);

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
