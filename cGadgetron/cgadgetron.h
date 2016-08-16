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
	void* cGT_getAcquisitionsData
		(void* ptr_acqs, unsigned int slice, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);

	void* cGT_reconstructImages(void* ptr_recon, void* ptr_input);
	void* cGT_reconstructedImages(void* ptr_recon);
	void* cGT_processImages(void* ptr_proc, void* ptr_input);
	void* cGT_selectImages(void* ptr_input, unsigned int inc, unsigned int off);
	void cGT_setImageToRealConversion(void* ptr_imgs, int type);
	void* cGT_imagesCopy(const void* ptr_imgs);
	void* cGT_writeImages
		(void* ptr_imgs, const char* out_file, const char* out_group);
	void* cGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num);
	void* cGT_imageTypes(const void* ptr_x);

	void cGT_getCoilDataDimensions(void* ptr_csms, int csm_num, PTR_INT ptr_dim);
	void cGT_getCoilData
		(void* ptr_csms, int csm_num, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
	void cGT_getCoilDataAbs(void* ptr_csms, int csm_num, PTR_DOUBLE ptr);
	void cGT_getImageDimensions(void* ptr_imgs, int img_num, PTR_INT ptr_dim);
	void cGT_getImageDataAsDoubleArray
		(void* ptr_imgs, int img_num, PTR_DOUBLE ptr_data);
	void cGT_getImageDataAsComplexArray
		(void* ptr_imgs, int img_num, PTR_DOUBLE ptr_data);

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
