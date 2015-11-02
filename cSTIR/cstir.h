#ifndef cSTIR_INTERFACE
#define cSTIR_INTERFACE

extern "C" {
	void* cSTIR_setParameter
		(void* ptr, const char* set, const char* name, const void* value);
	void* cSTIR_parameter(const void* ptr, const char* set, const char* name);
	void* cSTIR_newObject(const char* name);
	void cSTIR_deleteObject(void* ptr, const char* name);
	void* cSTIR_setupObject(const char* obj, void* ptr_obj);
	void* cSTIR_newReconstruction(const char* method, const char* filename);
	void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_reconstruct(void* ptr_r, void* ptr_i);
	void* cSTIR_update(void* ptr_r, void* ptr_i);
#ifndef CSTIR_FOR_MATLAB
	void cSTIR_getImageDimensions(const void* ptr, size_t pd);
	void cSTIR_getImageData(const void* ptr, size_t pd);
#else
	void cSTIR_getImageDimensions(const void* ptr, int* pd);
	void cSTIR_getImageData(const void* ptr, double* pd);
#endif
	void* cSTIR_voxels3DF(int nx, int ny, int nz,
		double sx, double sy, double sz, double x, double y, double z);
	void* cSTIR_imageFromVoxels(void* ptr_v);
	void* cSTIR_imageFromFile(const char* filename);
	void* cSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v);
	void* cSTIR_imagesDifference(void* first, void* second, int rimsize);

	void* newDataHandle();
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
}

#endif