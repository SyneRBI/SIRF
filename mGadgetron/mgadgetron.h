#ifndef CGADGETRON_TO_MATLAB_INTERFACE
#define CGADGETRON_TO_MATLAB_INTERFACE

#define CGADGETRON_FOR_MATLAB
#include "shrhelp.h"

#ifndef CGADGETRON_FOR_MATLAB
extern "C" {
#endif
EXPORTED_FUNCTION void* mGT_newObject(const char* name);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModel(const void* ptr_acqs, const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelForward(void* ptr_am, const void* ptr_imgs);
EXPORTED_FUNCTION 	void*	mGT_AcquisitionModelBackward(void* ptr_am, const void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelFwd(void* ptr_am, const void* ptr_imgs, void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_AcquisitionModelBwd(void* ptr_am, const void* ptr_imgs, void* ptr_acqs);
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFromFile(const char* file);
EXPORTED_FUNCTION 	void* mGT_ISMRMRDAcquisitionsFile(const char* file);
EXPORTED_FUNCTION 	void* mGT_newAcquisitionsContainer(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_acquisitionsNorm(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_acquisitionsDot(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_acquisitionsAxpby(double a, const void* ptr_x, double b, const void* ptr_y); 
EXPORTED_FUNCTION 	void* mGT_newImagesContainer(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_imagesCopy(const void* ptr_imgs);
EXPORTED_FUNCTION 	void* mGT_imagesNorm(const void* ptr_x);
EXPORTED_FUNCTION 	void* mGT_imagesDot(const void* ptr_x, const void* ptr_y);
EXPORTED_FUNCTION 	void* mGT_imagesAxpby(double a, const void* ptr_x, double b, const void* ptr_y); 
EXPORTED_FUNCTION 	void* mGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms);
EXPORTED_FUNCTION 	void* mGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION 	void* mGT_setGadgetProperty(void* ptr_g, const char* prop, const char* value);
EXPORTED_FUNCTION 	void* mGT_configGadgetChain(void* ptr_con, void* ptr_gc);
EXPORTED_FUNCTION 	void* mGT_registerHDFReceiver(void* ptr_con, const char* file, const char* group);
EXPORTED_FUNCTION 	void* mGT_registerImagesReceiver(void* ptr_con, void* ptr_img);
EXPORTED_FUNCTION 	void* mGT_reconstructImages(void* ptr_recon, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_reconstructedImages(void* ptr_recon);
EXPORTED_FUNCTION 	void* mGT_processImages(void* ptr_proc, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_processAcquisitions(void* ptr_proc, void* ptr_input);
EXPORTED_FUNCTION 	void* mGT_writeImages(void* ptr_imgs, const char* out_file, const char* out_group);
EXPORTED_FUNCTION 	int mGT_numImages(void* ptr_imgs);
#ifndef CGADGETRON_FOR_MATLAB
EXPORTED_FUNCTION 	void mGT_getImageDimensions(void* ptr_imgs, int im_num, size_t ptr_dim);
EXPORTED_FUNCTION 	void mGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, size_t ptr_data);
#else
EXPORTED_FUNCTION 	void mGT_getImageDimensions(void* ptr_imgs, int im_num, int* dim);
EXPORTED_FUNCTION 	void mGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, double* data);
#endif
EXPORTED_FUNCTION 	void* mGT_connect(void* ptr_con, const char* host, const char* port);
EXPORTED_FUNCTION 	void* mGT_sendConfigScript(void* ptr_con, const char* config);
EXPORTED_FUNCTION 	void* mGT_sendConfigFile(void* ptr_con, const char* file);
EXPORTED_FUNCTION 	void* mGT_sendParameters(void* ptr_con, const void* par);
EXPORTED_FUNCTION 	void* mGT_sendParametersString(void* ptr_con, const char* par);
EXPORTED_FUNCTION 	void* mGT_sendAcquisitions(void* ptr_con, void* ptr_dat);
EXPORTED_FUNCTION 	void* mGT_sendImages(void* ptr_con, void* ptr_img);
EXPORTED_FUNCTION 	void* mGT_disconnect(void* ptr_con);
EXPORTED_FUNCTION 	double mDoubleDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	double mDoubleReDataFromHandle(const void* ptr);
EXPORTED_FUNCTION 	double mDoubleImDataFromHandle(const void* ptr);
#ifndef CGADGETRON_FOR_MATLAB
}
#endif


#endif
