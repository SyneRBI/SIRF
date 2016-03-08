#ifndef GADGETRON_C_INTERFACE
#define GADGETRON_C_INTERFACE

extern "C" {
	void* cGT_newObject(const char* name);
	void* cGT_AcquisitionModel(const void* ptr_acqs);
	void* cGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
	void* cGT_AcquisitionModelFwd
		(void* ptr_am, const void* ptr_imgs, void* ptr_acqs);
	void* cGT_AcquisitionModelBwd
		(void* ptr_am, const void* ptr_imgs, void* ptr_acqs);
	void* cGT_ISMRMRDAcquisitionsFromFile(const char* file);
	void* cGT_ISMRMRDAcquisitionsFile(const char* file);
	void* cGT_newAcquisitionsContainer(const void* ptr_x);
	void* cGT_acquisitionsNorm(const void* ptr_x);
	void* cGT_acquisitionsDot(const void* ptr_x, const void* ptr_y);
	void* cGT_acquisitionsAxpby
		(double a, const void* ptr_x, double b, const void* ptr_y, void* ptr_z);
	void* cGT_imagesCopy(const void* ptr_imgs);
	void* cGT_imagesDot(const void* ptr_x, const void* ptr_y);
	void* cGT_acquisitionsProcessor();
	void* cGT_ISMRMRDatasetFromFile(const char* file, const char* group);
	void* cGT_readISMRMRDatasetHeader(void* ptr_data, void* ptr_head);
	void* cGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms);
	void* cGT_connect(void* ptr_con, const char* host, const char* port);
	void* cGT_sendConfigScript(void* ptr_con, const char* config);
	void* cGT_sendConfigFile(void* ptr_con, const char* file);
	void* cGT_sendParameters(void* ptr_con, const void* par);
	void* cGT_sendParametersString(void* ptr_con, const char* par);
	void* cGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
	void* cGT_setEndGadget(void* ptr_gc, const void* ptr_g);
	void* cGT_setGadgetProperty(void* ptr_g, const char* prop, const char* value);
	void* cGT_configGadgetChain(void* ptr_con, void* ptr_gc);
	void* cGT_registerHDFReceiver
		(void* ptr_con, const char* file, const char* group);
	void* cGT_registerImagesReceiver(void* ptr_con, void* ptr_img);
	void* cGT_reconstructImages(void* ptr_recon, void* ptr_input);
	void* cGT_reconstructedImagesList(void* ptr_recon);
	void* cGT_processImages(void* ptr_proc, void* ptr_input);
	void* cGT_processAcquisitions(void* ptr_proc, void* ptr_input);
	void* cGT_writeImages
		(void* ptr_imgs, const char* out_file, const char* out_group);
	int cGT_numImages(void* ptr_imgs);
#ifndef CGADGETRON_FOR_MATLAB
	void cGT_getImageDimensions(void* ptr_imgs, int im_num, size_t ptr_dim);
	void cGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, size_t ptr_data);
#else
	void cGT_getImageDimensions(void* ptr_imgs, int im_num, int* dim);
	void cGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, double* data);
#endif
	void* cGT_sendAcquisitions(void* ptr_con, void* ptr_dat);
	void* cGT_sendImages(void* ptr_con, void* ptr_img);
	void* cGT_disconnect(void* ptr_con);

	double doubleDataFromHandle(const void* ptr);
	double doubleReDataFromHandle(const void* ptr);
	double doubleImDataFromHandle(const void* ptr);
}

#endif
