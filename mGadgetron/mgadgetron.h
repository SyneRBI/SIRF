#ifndef CGADGETRON_TO_MATLAB_INTERFACE
#define CGADGETRON_TO_MATLAB_INTERFACE

#define CGADGETRON_FOR_MATLAB
#include "shrhelp.h"

EXPORTED_FUNCTION void* mGT_ISMRMRDatasetFromFile(const char* file, const char* group);
EXPORTED_FUNCTION void* mGT_readISMRMRDatasetHeader(void* ptr_data, void* ptr_head);
EXPORTED_FUNCTION void* mGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms);
EXPORTED_FUNCTION void* mGT_connect(void* ptr_con, const char* host, const char* port);
EXPORTED_FUNCTION void* mGT_sendConfigScript(void* ptr_con, const char* config);
EXPORTED_FUNCTION void* mGT_sendConfigFile(void* ptr_con, const char* file);
EXPORTED_FUNCTION void* mGT_sendParameters(void* ptr_con, const void* par);
EXPORTED_FUNCTION void* mGT_sendParametersString(void* ptr_con, const char* par);
EXPORTED_FUNCTION void* mGT_addReader(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION void* mGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION void* mGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r);
EXPORTED_FUNCTION void* mGT_configGadgetChain(void* ptr_con, void* ptr_gc);
EXPORTED_FUNCTION void* mGT_registerHDFReceiver(void* ptr_con, const char* file, const char* group);
EXPORTED_FUNCTION void* mGT_registerImagesReceiver(void* ptr_con, void* ptr_img);
EXPORTED_FUNCTION void* mGT_runMRIReconstruction(void* ptr_recon, void* ptr_input);
EXPORTED_FUNCTION void* mGT_reconstructedImagesList(void* ptr_recon);
EXPORTED_FUNCTION void* mGT_writeImages(void* ptr_imgs, const char* out_file, const char* out_group);
EXPORTED_FUNCTION int mGT_numImages(void* ptr_imgs);
#ifndef CGADGETRON_FOR_MATLAB
EXPORTED_FUNCTION void mGT_getImageDimensions(void* ptr_imgs, int im_num, size_t ptr_dim);
EXPORTED_FUNCTION void mGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, size_t ptr_data);
#else
EXPORTED_FUNCTION void mGT_getImageDimensions(void* ptr_imgs, int im_num, int* dim);
EXPORTED_FUNCTION void mGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, double* data);
#endif
EXPORTED_FUNCTION void* mGT_sendAcquisitions(void* ptr_con, void* ptr_dat);
EXPORTED_FUNCTION void* mGT_sendImages(void* ptr_con, void* ptr_img);
EXPORTED_FUNCTION void* mGT_disconnect(void* ptr_con);
EXPORTED_FUNCTION void* mNewObject(const char* name);
EXPORTED_FUNCTION void* mCopyOfObject(void* ptr);
EXPORTED_FUNCTION void mDeleteObject(void* ptr);
EXPORTED_FUNCTION void* mNewDataHandle();
EXPORTED_FUNCTION void mDeleteDataHandle(void* ptr) ;
EXPORTED_FUNCTION int mExecutionStatus(const void* ptr);
EXPORTED_FUNCTION const char* mExecutionError(const void* ptr);
EXPORTED_FUNCTION const char* mExecutionErrorFile(const void* ptr);
EXPORTED_FUNCTION int mExecutionErrorLine(const void* ptr);

//EXPORTED_FUNCTION void* mNewMexPrinter();
//EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr);
#endif
