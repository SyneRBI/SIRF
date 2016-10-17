#ifndef CSTIR_TO_MATLAB_INTERFACE
#define CSTIR_TO_MATLAB_INTERFACE

#define CSTIR_FOR_MATLAB
#include "shrhelp.h"

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
EXPORTED_FUNCTION  void* mSTIR_newObject(const char* name);
EXPORTED_FUNCTION 	void* mSTIR_objectFromFile(const char* name, const char* filename);
EXPORTED_FUNCTION 	void* mSTIR_setParameter (void* ptr, const char* obj, const char* name, const void* value);
EXPORTED_FUNCTION 	void* mSTIR_parameter(const void* ptr, const char* obj, const char* name);
EXPORTED_FUNCTION 	void* mSTIR_setupObject(const char* obj, void* ptr_obj);
EXPORTED_FUNCTION 	void* mSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelSetup(void* ptr_am, const char* templ, void* ptr_im);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelForward (void* ptr_am, const char* datafile, void* ptr_dt, void* ptr_im);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelBackward(void* ptr_am, void* ptr_ad, void* ptr_im);
EXPORTED_FUNCTION 	void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_value(void* ptr_f, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_gradient(void* ptr_f, void* ptr_i, int subset);
EXPORTED_FUNCTION 	void mSTIR_getImageDimensions(const void* ptr, PTR_INT pd);
EXPORTED_FUNCTION 	void mSTIR_getImageData(const void* ptr, PTR_DOUBLE pd);
EXPORTED_FUNCTION 	void mSTIR_setImageData(const void* ptr_im, size_t ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_voxels3DF(int nx, int ny, int nz, double sx, double sy, double sz, double x, double y, double z);
EXPORTED_FUNCTION 	void* mSTIR_imageFromVoxels(void* ptr_v);
EXPORTED_FUNCTION 	void* mSTIR_imageFromImage(void* ptr_v);
EXPORTED_FUNCTION 	void mSTIR_fillImage(void* ptr_i, double v);
EXPORTED_FUNCTION 	void* mSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v);
EXPORTED_FUNCTION 	void* mSTIR_imagesDifference(void* first, void* second, int rimsize);
EXPORTED_FUNCTION 	void* mNewTextPrinter(const char* stream);
EXPORTED_FUNCTION 	void* mNewTextWriter(const char* stream);
EXPORTED_FUNCTION 	void mOpenChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION 	void mCloseChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION 	void mSetWriter(void* ptr_w, int channel);
EXPORTED_FUNCTION 	void mResetWriter();
EXPORTED_FUNCTION 	void mDeleteTextPrinter(void* ptr);
EXPORTED_FUNCTION 	void mDeleteTextWriter(void* ptr_w);
#ifndef CSTIR_FOR_MATLAB
}
#endif
EXPORTED_FUNCTION void* mNewMexPrinter();
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr);

#endif
