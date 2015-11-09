#ifndef CSTIR_INTERFACE
#define CSTIR_INTERFACE

void* cSTIR_newObject(const char* name);
void cSTIR_deleteObject(void* ptr, const char* name);
void* cSTIR_copyOfObject(const char* name, void* ptr);
void* cSTIR_setParameter(void* ptr, const char* set, const char* name, const void* value);
void* cSTIR_parameter(const void* ptr, const char* set, const char* name);
void* cSTIR_setupObject(const char* obj, void* ptr_obj);
void* cSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d);
void* cSTIR_newReconstruction(const char* method, const char* filename);
void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
void* cSTIR_runReconstruction(void* ptr_r, void* ptr_i);
void* cSTIR_updateReconstruction(void* ptr_r, void* ptr_i);
void cSTIR_getImageDimensions(const void* ptr, size_t pd);
void cSTIR_getImageData(const void* ptr, size_t pd);
void cSTIR_getImageDimensions(const void* ptr, int* pd);
void cSTIR_getImageData(const void* ptr, double* pd);
void* cSTIR_voxels3DF(int nx, int ny, int nz,double sx, double sy, double sz, double x, double y, double z);
void* cSTIR_imageFromVoxels(void* ptr_v);
void* cSTIR_imageFromImage(void* ptr_v);
void* cSTIR_imageFromFile(const char* filename);
void cSTIR_fillImage(void* ptr_i, double v);
void* cSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v);
void* cSTIR_imagesDifference(void* first, void* second, int rimsize);
void* newDataHandle();
void* refDataHandle(void* ptr);
void* charDataHandle(const char* s);
void* intDataHandle(int i);
void* floatDataHandle(float i);
void* doubleDataHandle(double i);
char* charDataFromHandle(const void* ptr);
int intDataFromHandle(const void* ptr);
float floatDataFromHandle(const void* ptr);
double doubleDataFromHandle(const void* ptr);
void deleteDataHandle(void* ptr);
int executionStatus(const void* ptr);
const char* executionError(const void* ptr);
const char* executionErrorFile(const void* ptr);
int executionErrorLine(const void* ptr);
void* newTextPrinter(const char* stream);
void* newTextWriter(const char* stream);
void openChannel(int channel, void* ptr_w);
void closeChannel(int channel);
void setWriter(void* ptr_w, int channel);
void resetWriter();
void deleteTextPrinter(void* ptr);
void deleteTextWriter(void* ptr_w);

#endif
