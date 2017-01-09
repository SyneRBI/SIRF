#ifndef CGADGETRON_TO_MATLAB_INTERFACE
#define CGADGETRON_TO_MATLAB_INTERFACE

#define CGADGETRON_FOR_MATLAB
#include "shrhelp.h"

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
EXPORTED_FUNCTION 	void* mGT_appendCSM (void* ptr_csms, int nx, int ny, int nz, int nc,  PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_setCSMs(void* ptr_am, const void* ptr_csms);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFromFile(const char* file);
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFile(const char* file);
EXPORTED_FUNCTION 	void* mGT_processAcquisitions(void* ptr_proc, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_acquisitionFromContainer(void* ptr_acqs, unsigned int acq_num);
EXPORTED_FUNCTION 	void* mGT_orderAcquisitions(void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsDimensions(void* ptr_acqs, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsFlags(void* ptr_acqs, unsigned int n, PTR_INT ptr_f);
EXPORTED_FUNCTION 	void* mGT_getAcquisitionsData (void* ptr_acqs, unsigned int slice, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
EXPORTED_FUNCTION 	void* mGT_setAcquisitionsData (void* ptr_acqs, unsigned int na, unsigned int nc, unsigned int ns, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
EXPORTED_FUNCTION 	void* mGT_reconstructImages(void* ptr_recon, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_reconstructedImages(void* ptr_recon);
EXPORTED_FUNCTION 	void* mGT_processImages(void* ptr_proc, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_selectImages(void* ptr_input, unsigned int inc, unsigned int off);
EXPORTED_FUNCTION 	void mGT_setImageToRealConversion(void* ptr_imgs, int type);
EXPORTED_FUNCTION 	void* mGT_imagesCopy(const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_writeImages (void* ptr_imgs, const char* out_file, const char* out_group);
EXPORTED_FUNCTION 	void* mGT_imageWrapFromContainer(void* ptr_imgs, unsigned int img_num);
EXPORTED_FUNCTION 	void* mGT_imageTypes(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_imageDataType(const void* ptr_x, int im_num);
EXPORTED_FUNCTION 	void mGT_getCoilDataDimensions(void* ptr_csms, int csm_num, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void mGT_getCoilData (void* ptr_csms, int csm_num, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
EXPORTED_FUNCTION 	void mGT_getCoilDataAbs(void* ptr_csms, int csm_num, PTR_DOUBLE ptr);
EXPORTED_FUNCTION 	void mGT_getImageDimensions(void* ptr_imgs, int img_num, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void mGT_getImageDataAsDoubleArray (void* ptr_imgs, int img_num, PTR_DOUBLE ptr_data);
EXPORTED_FUNCTION 	void mGT_getImageDataAsComplexArray (void* ptr_imgs, int img_num, PTR_DOUBLE ptr_data);
EXPORTED_FUNCTION 	void mGT_getImagesDataAsDoubleArray(void* ptr_imgs, PTR_DOUBLE ptr_data);
EXPORTED_FUNCTION 	void mGT_getImagesDataAsComplexArray (void* ptr_imgs, PTR_DOUBLE ptr_re, PTR_DOUBLE ptr_im);
EXPORTED_FUNCTION 	void* mGT_dataItems(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_norm(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_dot(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_axpby( double ar, double ai, const void* ptr_x, double br, double bi, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_setGadgetProperty(void* ptr_g, const char* prop, const char* value);
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
