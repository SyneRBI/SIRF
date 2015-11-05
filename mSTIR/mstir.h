#ifndef CSTIR_TO_MATLAB_INTERFACE
#define CSTIR_TO_MATLAB_INTERFACE

#define CSTIR_FOR_MATLAB
#include "shrhelp.h"

EXPORTED_FUNCTION void* mSTIR_setParameter(void* ptr, const char* set, const char* name, const void* value);
EXPORTED_FUNCTION void* mSTIR_parameter(const void* ptr, const char* set, const char* name);
EXPORTED_FUNCTION void* mSTIR_newObject(const char* name);
EXPORTED_FUNCTION void mSTIR_deleteObject(void* ptr, const char* name);
EXPORTED_FUNCTION void* mSTIR_setupObject(const char* obj, void* ptr_obj);
EXPORTED_FUNCTION void* mSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d);
EXPORTED_FUNCTION void* mSTIR_newReconstruction(const char* method, const char* filename);
EXPORTED_FUNCTION void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i);
#ifndef CSTIR_FOR_MATLAB
EXPORTED_FUNCTION void mSTIR_getImageDimensions(const void* ptr, size_t pd);
EXPORTED_FUNCTION void mSTIR_getImageData(const void* ptr, size_t pd);
#else
EXPORTED_FUNCTION void mSTIR_getImageDimensions(const void* ptr, int* pd);
EXPORTED_FUNCTION void mSTIR_getImageData(const void* ptr, double* pd);
#endif
EXPORTED_FUNCTION void* mSTIR_voxels3DF(int nx, int ny, int nz,double sx, double sy, double sz, double x, double y, double z);
EXPORTED_FUNCTION void* mSTIR_imageFromVoxels(void* ptr_v);
EXPORTED_FUNCTION void* mSTIR_imageFromImage(void* ptr_v);
EXPORTED_FUNCTION void* mSTIR_imageFromFile(const char* filename);
EXPORTED_FUNCTION void mSTIR_fillImage(void* ptr_i, double v);
EXPORTED_FUNCTION void* mSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v);
EXPORTED_FUNCTION void* mSTIR_imagesDifference(void* first, void* second, int rimsize);
EXPORTED_FUNCTION void* mNewDataHandle();
EXPORTED_FUNCTION void* mCharDataHandle(const char* s);
EXPORTED_FUNCTION void* mIntDataHandle(int i);
EXPORTED_FUNCTION void* mFloatDataHandle(float i);
EXPORTED_FUNCTION void* mDoubleDataHandle(double i);
EXPORTED_FUNCTION char* mCharDataFromHandle(const void* ptr);
EXPORTED_FUNCTION int mIntDataFromHandle(const void* ptr);
EXPORTED_FUNCTION float mFloatDataFromHandle(const void* ptr);
EXPORTED_FUNCTION double mDoubleDataFromHandle(const void* ptr);
EXPORTED_FUNCTION void mDeleteDataHandle(void* ptr);
EXPORTED_FUNCTION int mExecutionStatus(const void* ptr);
EXPORTED_FUNCTION const char* mExecutionError(const void* ptr);
EXPORTED_FUNCTION const char* mExecutionErrorFile(const void* ptr);
EXPORTED_FUNCTION int mExecutionErrorLine(const void* ptr);
EXPORTED_FUNCTION void* mNewTextPrinter(const char* stream);
EXPORTED_FUNCTION void* mNewTextWriter(const char* stream);
EXPORTED_FUNCTION void mOpenChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION void mCloseChannel(int channel);
EXPORTED_FUNCTION void mSetWriter(void* ptr_w, int channel);
EXPORTED_FUNCTION void mResetWriter();
EXPORTED_FUNCTION void mDeleteTextPrinter(void* ptr);
EXPORTED_FUNCTION void mDeleteTextWriter(void* ptr_w);

EXPORTED_FUNCTION void* mNewMexPrinter();
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr);
#endif
