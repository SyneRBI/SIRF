#ifndef cSTIR_INTERFACE
#define cSTIR_INTERFACE

extern "C" {

	// Common STIR Object methods
	void* cSTIR_newObject(const char* name);
	void cSTIR_deleteObject(void* ptr);
	void* cSTIR_objectFromFile(const char* name, const char* filename);
	void* cSTIR_copyOfObject(void* ptr);
	void* cSTIR_setParameter
		(void* ptr, const char* obj, const char* name, const void* value);
	void* cSTIR_parameter(const void* ptr, const char* obj, const char* name);
	void* cSTIR_setupObject(const char* obj, void* ptr_obj);

	// DataProcessor methods
	void* cSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d);

	// Acquisition model methods
	void* cSTIR_acquisitionModelSetup(void* ptr_am, const char* templ, void* ptr_im);
	void* cSTIR_acquisitionModelForward
		(void* ptr_am, const char* datafile, void* ptr_dt, void* ptr_im);
	void* cSTIR_acquisitionModelBackward(void* ptr_am, void* ptr_ad, void* ptr_im);

	// Reconstruction methods
	void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_runReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_updateReconstruction(void* ptr_r, void* ptr_i);

	// Image methods
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
	void* cSTIR_imageFromImage(void* ptr_v);
	void cSTIR_fillImage(void* ptr_i, double v);
	void* cSTIR_addShape(void* ptr_i, void* ptr_v, void* ptr_s, float v);
	void* cSTIR_imagesDifference(void* first, void* second, int rimsize);

	// DataHandle methods
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

	// ExecutionStatus methods
	int executionStatus(const void* ptr);
	const char* executionError(const void* ptr);
	const char* executionErrorFile(const void* ptr);
	int executionErrorLine(const void* ptr);

	// TextWriter methods
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