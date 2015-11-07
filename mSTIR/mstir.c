#include <mex.h>
#define EXPORT_FCNS
#define CSTIR_FOR_MATLAB
#include "matrix.h"
#include "shrhelp.h"
#include "cstir.h"

EXPORTED_FUNCTION void* mSTIR_newObject(const char* name) {
	return cSTIR_newObject(name);
}
EXPORTED_FUNCTION void mSTIR_deleteObject(void* ptr, const char* name) {
	cSTIR_deleteObject(ptr, name);
}
EXPORTED_FUNCTION void* mSTIR_setParameter(void* ptr, const char* set, const char* name, const void* value) {
	return cSTIR_setParameter(ptr, set, name, value);
}
EXPORTED_FUNCTION void* mSTIR_parameter(const void* ptr, const char* set, const char* name) {
	return cSTIR_parameter(ptr, set, name);
}
EXPORTED_FUNCTION void* mSTIR_setupObject(const char* obj, void* ptr_obj) {
	return cSTIR_setupObject(obj, ptr_obj);
}
EXPORTED_FUNCTION void* mSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d) {
	return cSTIR_applyDataProcessor(ptr_p, ptr_d);
}
EXPORTED_FUNCTION void* mSTIR_newReconstruction(const char* method, const char* filename) {
	return cSTIR_newReconstruction(method, filename);
}
EXPORTED_FUNCTION void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_setupReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_runReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_updateReconstruction(ptr_r, ptr_i);
}
#ifndef CSTIR_FOR_MATLAB
EXPORTED_FUNCTION void mSTIR_getImageDimensions(const void* ptr, size_t pd) {
	cSTIR_getImageDimensions(ptr, pd);
}
EXPORTED_FUNCTION void mSTIR_getImageData(const void* ptr, size_t pd) {
	cSTIR_getImageData(ptr, pd);
}
#else
EXPORTED_FUNCTION void mSTIR_getImageDimensions(const void* ptr, int* pd) {
	cSTIR_getImageDimensions(ptr, pd);
}
EXPORTED_FUNCTION void mSTIR_getImageData(const void* ptr, double* pd) {
	cSTIR_getImageData(ptr, pd);
}
#endif
EXPORTED_FUNCTION void* mSTIR_voxels3DF(int nx, int ny, int nz,double sx, double sy, double sz, double x, double y, double z) {
	return cSTIR_voxels3DF(nx, ny, nz, sx, sy, sz, x, y, z);
}
EXPORTED_FUNCTION void* mSTIR_imageFromVoxels(void* ptr_v) {
	return cSTIR_imageFromVoxels(ptr_v);
}
EXPORTED_FUNCTION void* mSTIR_imageFromImage(void* ptr_v) {
	return cSTIR_imageFromImage(ptr_v);
}
EXPORTED_FUNCTION void* mSTIR_imageFromFile(const char* filename) {
	return cSTIR_imageFromFile(filename);
}
EXPORTED_FUNCTION void mSTIR_fillImage(void* ptr_i, double v) {
	cSTIR_fillImage(ptr_i, v);
}
EXPORTED_FUNCTION void* mSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v) {
	return cSTIR_addShape(ptr_i, ptr_v, ptr_s, v);
}
EXPORTED_FUNCTION void* mSTIR_imagesDifference(void* first, void* second, int rimsize) {
	return cSTIR_imagesDifference(first, second, rimsize);
}
EXPORTED_FUNCTION void* mNewDataHandle() {
	return newDataHandle();
}
EXPORTED_FUNCTION void* mRefDataHandle(void* ptr) {
	return refDataHandle(ptr);
}
EXPORTED_FUNCTION void* mCharDataHandle(const char* s) {
	return charDataHandle(s);
}
EXPORTED_FUNCTION void* mIntDataHandle(int i) {
	return intDataHandle(i);
}
EXPORTED_FUNCTION void* mFloatDataHandle(float i) {
	return floatDataHandle(i);
}
EXPORTED_FUNCTION void* mDoubleDataHandle(double i) {
	return doubleDataHandle(i);
}
EXPORTED_FUNCTION char* mCharDataFromHandle(const void* ptr) {
	return charDataFromHandle(ptr);
}
EXPORTED_FUNCTION int mIntDataFromHandle(const void* ptr) {
	return intDataFromHandle(ptr);
}
EXPORTED_FUNCTION float mFloatDataFromHandle(const void* ptr) {
	return floatDataFromHandle(ptr);
}
EXPORTED_FUNCTION double mDoubleDataFromHandle(const void* ptr) {
	return doubleDataFromHandle(ptr);
}
EXPORTED_FUNCTION void mDeleteDataHandle(void* ptr) {
	deleteDataHandle(ptr);
}
EXPORTED_FUNCTION int mExecutionStatus(const void* ptr) {
	return executionStatus(ptr);
}
EXPORTED_FUNCTION const char* mExecutionError(const void* ptr) {
	return executionError(ptr);
}
EXPORTED_FUNCTION const char* mExecutionErrorFile(const void* ptr) {
	return executionErrorFile(ptr);
}
EXPORTED_FUNCTION int mExecutionErrorLine(const void* ptr) {
	return executionErrorLine(ptr);
}
EXPORTED_FUNCTION void* mNewTextPrinter(const char* stream) {
	return newTextPrinter(stream);
}
EXPORTED_FUNCTION void* mNewTextWriter(const char* stream) {
	return newTextWriter(stream);
}
EXPORTED_FUNCTION void mOpenChannel(int channel, void* ptr_w) {
	openChannel(channel, ptr_w);
}
EXPORTED_FUNCTION void mCloseChannel(int channel) {
	closeChannel(channel);
}
EXPORTED_FUNCTION void mSetWriter(void* ptr_w, int channel) {
	setWriter(ptr_w, channel);
}
EXPORTED_FUNCTION void mResetWriter() {
	resetWriter();
}
EXPORTED_FUNCTION void mDeleteTextPrinter(void* ptr) {
	deleteTextPrinter(ptr);
}
EXPORTED_FUNCTION void mDeleteTextWriter(void* ptr_w) {
	deleteTextWriter(ptr_w);
}

EXPORTED_FUNCTION void* mNewMexPrinter() {
  return newMexPrinter();
}
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr) {
  deleteMexPrinter(ptr);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
