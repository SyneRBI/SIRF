#include <mex.h>
#define EXPORT_FCNS
#define CSTIR_FOR_MATLAB
#include "matrix.h"
#include "shrhelp.h"
#include "cstir.h"

#ifndef CSTIR_FOR_MATLAB
#define PTR_INT size_t
#define PTR_FLOAT size_t
#define PTR_DOUBLE size_t
extern "C" {
#else

#define PTR_INT int*

#define PTR_FLOAT float*

#define PTR_DOUBLE double*

#endif
EXPORTED_FUNCTION  void* mSTIR_newObject(const char* name) {
	return cSTIR_newObject(name);
}
EXPORTED_FUNCTION 	void mSTIR_deleteObject(void* ptr) {
	cSTIR_deleteObject(ptr);
}
EXPORTED_FUNCTION 	void* mSTIR_objectFromFile(const char* name, const char* filename) {
	return cSTIR_objectFromFile(name, filename);
}
EXPORTED_FUNCTION 	void* mSTIR_copyOfObject(void* ptr) {
	return cSTIR_copyOfObject(ptr);
}
EXPORTED_FUNCTION 	void* mSTIR_setParameter (void* ptr, const char* obj, const char* name, const void* value) {
	return cSTIR_setParameter (ptr, obj, name, value);
}
EXPORTED_FUNCTION 	void* mSTIR_parameter(const void* ptr, const char* obj, const char* name) {
	return cSTIR_parameter(ptr, obj, name);
}
EXPORTED_FUNCTION 	void* mSTIR_setupObject(const char* obj, void* ptr_obj) {
	return cSTIR_setupObject(obj, ptr_obj);
}
EXPORTED_FUNCTION 	void* mSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d) {
	return cSTIR_applyDataProcessor(ptr_p, ptr_d);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelSetup(void* ptr_am, const char* templ, void* ptr_im) {
	return cSTIR_acquisitionModelSetup(ptr_am, templ, ptr_im);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelForward (void* ptr_am, const char* datafile, void* ptr_dt, void* ptr_im) {
	return cSTIR_acquisitionModelForward (ptr_am, datafile, ptr_dt, ptr_im);
}
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelBackward(void* ptr_am, void* ptr_ad, void* ptr_im) {
	return cSTIR_acquisitionModelBackward(ptr_am, ptr_ad, ptr_im);
}
EXPORTED_FUNCTION 	void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_setupReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_runReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i) {
	return cSTIR_updateReconstruction(ptr_r, ptr_i);
}
EXPORTED_FUNCTION 	void mSTIR_getImageDimensions(const void* ptr, PTR_INT pd) {
	cSTIR_getImageDimensions(ptr, pd);
}
EXPORTED_FUNCTION 	void mSTIR_getImageData(const void* ptr, PTR_DOUBLE pd) {
	cSTIR_getImageData(ptr, pd);
}
EXPORTED_FUNCTION 	void* mSTIR_voxels3DF(int nx, int ny, int nz, double sx, double sy, double sz, double x, double y, double z) {
	return cSTIR_voxels3DF(nx, ny, nz, sx, sy, sz, x, y, z);
}
EXPORTED_FUNCTION 	void* mSTIR_imageFromVoxels(void* ptr_v) {
	return cSTIR_imageFromVoxels(ptr_v);
}
EXPORTED_FUNCTION 	void* mSTIR_imageFromImage(void* ptr_v) {
	return cSTIR_imageFromImage(ptr_v);
}
EXPORTED_FUNCTION 	void mSTIR_fillImage(void* ptr_i, double v) {
	cSTIR_fillImage(ptr_i, v);
}
EXPORTED_FUNCTION 	void* mSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v) {
	return cSTIR_addShape(ptr_i, ptr_v, ptr_s, v);
}
EXPORTED_FUNCTION 	void* mSTIR_imagesDifference(void* first, void* second, int rimsize) {
	return cSTIR_imagesDifference(first, second, rimsize);
}
EXPORTED_FUNCTION 	void* mNewTextPrinter(const char* stream) {
	return newTextPrinter(stream);
}
EXPORTED_FUNCTION 	void* mNewTextWriter(const char* stream) {
	return newTextWriter(stream);
}
EXPORTED_FUNCTION 	void mOpenChannel(int channel, void* ptr_w) {
	openChannel(channel, ptr_w);
}
EXPORTED_FUNCTION 	void mCloseChannel(int channel, void* ptr_w) {
	closeChannel(channel, ptr_w);
}
EXPORTED_FUNCTION 	void mSetWriter(void* ptr_w, int channel) {
	setWriter(ptr_w, channel);
}
EXPORTED_FUNCTION 	void mResetWriter() {
	resetWriter();
}
EXPORTED_FUNCTION 	void mDeleteTextPrinter(void* ptr) {
	deleteTextPrinter(ptr);
}
EXPORTED_FUNCTION 	void mDeleteTextWriter(void* ptr_w) {
	deleteTextWriter(ptr_w);
}
#ifndef CSTIR_FOR_MATLAB
}
#endif

void* newMexPrinter();
void deleteMexPrinter(void* ptr);
EXPORTED_FUNCTION void* mNewMexPrinter() {
  return newMexPrinter();
}
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr) {
  deleteMexPrinter(ptr);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {}
