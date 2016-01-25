#ifndef GADGETRON_C_INTERFACE
#define GADGETRON_C_INTERFACE

extern "C" {
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
	void* cGT_configGadgetChain(void* ptr_con, void* ptr_gc);
	void* cGT_registerHDFReceiver
		(void* ptr_con, const char* file, const char* group);
	void* cGT_registerImagesReceiver(void* ptr_con, void* ptr_img);
	void* cGT_runMRIReconstruction(void* ptr_recon, void* ptr_input);
	void* cGT_reconstructedImagesList(void* ptr_recon);
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

	void* newObject(const char* name);
	void* copyOfObject(void* ptr);
	void deleteObject(void* ptr);

	void* newDataHandle();
	void deleteDataHandle(void* ptr) ;
	int executionStatus(const void* ptr);
	const char* executionError(const void* ptr);
	const char* executionErrorFile(const void* ptr);
	int executionErrorLine(const void* ptr);
}

#endif
