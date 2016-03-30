#ifndef GADGETRON_C_INTERFACE
#define GADGETRON_C_INTERFACE

#ifndef CGADGETRON_FOR_MATLAB
extern "C" {
#endif
	void* cGT_newObject(const char* name);
	void* cGT_parameter(void* ptr, const char* obj, const char* name);

	void* cGT_CoilSensitivitiesFromFile(const char* file);

	void* cGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs);
	void* cGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
	void* cGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs);

	void* cGT_ISMRMRDAcquisitionsFromFile(const char* file);
	void* cGT_ISMRMRDAcquisitionsFile(const char* file);
	void* cGT_processAcquisitions(void* ptr_proc, void* ptr_input);
	void* cGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num);

	void* cGT_reconstructImages(void* ptr_recon, void* ptr_input);
	void* cGT_reconstructedImages(void* ptr_recon);
	void* cGT_processImages(void* ptr_proc, void* ptr_input);

	void* cGT_imagesCopy(const void* ptr_imgs);
	void* cGT_writeImages
		(void* ptr_imgs, const char* out_file, const char* out_group);
	void* cGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num);

#ifndef CGADGETRON_FOR_MATLAB
	void cGT_getCSMDimensions(void* ptr_csms, int csm_num, size_t ptr_dim);
	void
		cGT_getCSMData(void* ptr_csms, int csm_num, size_t ptr_re, size_t ptr_im);
	void cGT_getCSMDataAbs(void* ptr_csms, int csm_num, size_t ptr);
	void cGT_getImageDimensions(void* ptr_imgs, int img_num, size_t ptr_dim);
	void cGT_getImageDataAsDoubleArray
		(void* ptr_imgs, int img_num, size_t ptr_data);
	void cGT_getImageDataAsComplexArray
		(void* ptr_imgs, int img_num, size_t ptr_data);
#else
	void cGT_getCSMDimensions(void* ptr_csms, int csm_num, int* dim);
	void
		cGT_getCSMData(void* ptr_csms, int csm_num, double* re, double* im);
	void cGT_getCSMDataAbs(void* ptr_csms, int csm_num, double* v);
	void cGT_getImageDimensions(void* ptr_imgs, int img_num, int* dim);
	void cGT_getImageDataAsDoubleArray(void* ptr_imgs, int img_num, double* data);
	void cGT_getImageDataAsCmplxArray
		(void* ptr_imgs, int img_num, double* re, double* im);
#endif

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
	void* cGT_configGadgetChain(void* ptr_con, void* ptr_gc);
	void* cGT_registerHDFReceiver
		(void* ptr_con, const char* file, const char* group);
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
