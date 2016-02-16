#include <mex.h>
#define EXPORT_FCNS
#define CGADGETRON_FOR_MATLAB
#include "matrix.h"
#include "shrhelp.h"
#include "cgadgetron.h"

EXPORTED_FUNCTION void* mGT_newObject(const char* name) {
	return cGT_newObject(name);
}
EXPORTED_FUNCTION void* mGT_ISMRMRDAcquisitionsFromFile(const char* file) {
	return cGT_ISMRMRDAcquisitionsFromFile(file);
}
EXPORTED_FUNCTION void* mGT_acquisitionsProcessor(const char* file) {
	return cGT_acquisitionsProcessor(file);
}
EXPORTED_FUNCTION void* mGT_ISMRMRDatasetFromFile(const char* file, const char* group) {
	return cGT_ISMRMRDatasetFromFile(file, group);
}
EXPORTED_FUNCTION void* mGT_readISMRMRDatasetHeader(void* ptr_data, void* ptr_head) {
	return cGT_readISMRMRDatasetHeader(ptr_data, ptr_head);
}
EXPORTED_FUNCTION void* mGT_setConnectionTimeout(void* ptr_con, unsigned int timeout_ms) {
	return cGT_setConnectionTimeout(ptr_con, timeout_ms);
}
EXPORTED_FUNCTION void* mGT_connect(void* ptr_con, const char* host, const char* port) {
	return cGT_connect(ptr_con, host, port);
}
EXPORTED_FUNCTION void* mGT_sendConfigScript(void* ptr_con, const char* config) {
	return cGT_sendConfigScript(ptr_con, config);
}
EXPORTED_FUNCTION void* mGT_sendConfigFile(void* ptr_con, const char* file) {
	return cGT_sendConfigFile(ptr_con, file);
}
EXPORTED_FUNCTION void* mGT_sendParameters(void* ptr_con, const void* par) {
	return cGT_sendParameters(ptr_con, par);
}
EXPORTED_FUNCTION void* mGT_sendParametersString(void* ptr_con, const char* par) {
	return cGT_sendParametersString(ptr_con, par);
}
EXPORTED_FUNCTION void* mGT_addReader(void* ptr_gc, const char* id, const void* ptr_r) {
	return cGT_addReader(ptr_gc, id, ptr_r);
}
EXPORTED_FUNCTION void* mGT_addWriter(void* ptr_gc, const char* id, const void* ptr_r) {
	return cGT_addWriter(ptr_gc, id, ptr_r);
}
EXPORTED_FUNCTION void* mGT_addGadget(void* ptr_gc, const char* id, const void* ptr_r) {
	return cGT_addGadget(ptr_gc, id, ptr_r);
}
EXPORTED_FUNCTION void* mGT_setEndGadget(void* ptr_gc, const void* ptr_g) {
	return cGT_setEndGadget(ptr_gc, ptr_g);
}
EXPORTED_FUNCTION void* mGT_setGadgetProperty(void* ptr_g, const char* prop, const char* value) {
	return cGT_setGadgetProperty(ptr_g, prop, value);
}
EXPORTED_FUNCTION void* mGT_configGadgetChain(void* ptr_con, void* ptr_gc) {
	return cGT_configGadgetChain(ptr_con, ptr_gc);
}
EXPORTED_FUNCTION void* mGT_registerHDFReceiver(void* ptr_con, const char* file, const char* group) {
	return cGT_registerHDFReceiver(ptr_con, file, group);
}
EXPORTED_FUNCTION void* mGT_registerImagesReceiver(void* ptr_con, void* ptr_img) {
	return cGT_registerImagesReceiver(ptr_con, ptr_img);
}
EXPORTED_FUNCTION void* mGT_reconstructImages(void* ptr_recon, void* ptr_input) {
	return cGT_reconstructImages(ptr_recon, ptr_input);
}
EXPORTED_FUNCTION void* mGT_reconstructedImagesList(void* ptr_recon) {
	return cGT_reconstructedImagesList(ptr_recon);
}
EXPORTED_FUNCTION void* mGT_processImages(void* ptr_proc, void* ptr_input) {
	return cGT_processImages(ptr_proc, ptr_input);
}
EXPORTED_FUNCTION void* mGT_processAcquisitions(void* ptr_proc, void* ptr_input) {
	return cGT_processAcquisitions(ptr_proc, ptr_input);
}
EXPORTED_FUNCTION void* mGT_writeImages(void* ptr_imgs, const char* out_file, const char* out_group) {
	return cGT_writeImages(ptr_imgs, out_file, out_group);
}
EXPORTED_FUNCTION int mGT_numImages(void* ptr_imgs) {
	return cGT_numImages(ptr_imgs);
}
#ifndef CGADGETRON_FOR_MATLAB
EXPORTED_FUNCTION void mGT_getImageDimensions(void* ptr_imgs, int im_num, size_t ptr_dim) {
	cGT_getImageDimensions(ptr_imgs, im_num, ptr_dim);
}
EXPORTED_FUNCTION void mGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, size_t ptr_data) {
	cGT_getImageDataAsDoubleArray(ptr_imgs, im_num, ptr_data);
}
#else
EXPORTED_FUNCTION void mGT_getImageDimensions(void* ptr_imgs, int im_num, int* dim) {
	cGT_getImageDimensions(ptr_imgs, im_num, dim);
}
EXPORTED_FUNCTION void mGT_getImageDataAsDoubleArray(void* ptr_imgs, int im_num, double* data) {
	cGT_getImageDataAsDoubleArray(ptr_imgs, im_num, data);
}
#endif
EXPORTED_FUNCTION void* mGT_sendAcquisitions(void* ptr_con, void* ptr_dat) {
	return cGT_sendAcquisitions(ptr_con, ptr_dat);
}
EXPORTED_FUNCTION void* mGT_sendImages(void* ptr_con, void* ptr_img) {
	return cGT_sendImages(ptr_con, ptr_img);
}
EXPORTED_FUNCTION void* mGT_disconnect(void* ptr_con) {
	return cGT_disconnect(ptr_con);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
