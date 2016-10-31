#ifndef cSTIR_INTERFACE
#define cSTIR_INTERFACE

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

	// Common STIR Object methods
	void* cSTIR_newObject(const char* name);
	void* cSTIR_objectFromFile(const char* name, const char* filename);
	void* cSTIR_setParameter
		(void* ptr, const char* obj, const char* name, const void* value);
	void* cSTIR_parameter(const void* ptr, const char* obj, const char* name);
	//void* cSTIR_setupObject(const char* obj, void* ptr_obj);

	// DataProcessor methods
	void* cSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d);

	// Acquisition model methods
	void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im);
	void* cSTIR_acquisitionModelFwd
		(void* ptr_am, void* ptr_im, const char* datafile);
	void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad);

	// Acquisition data methods
	void* cSTIR_acquisitionsDataFromTemplate(void* ptr_t);
	void* cSTIR_getAcquisitionsDimensions(const void* ptr_acq, PTR_INT ptr_dim);
	void* cSTIR_getAcquisitionsData(const void* ptr_acq, PTR_DOUBLE ptr_data);
 	void cSTIR_setAcquisitionsData(void* ptr_acq, PTR_DOUBLE ptr_data);
	void cSTIR_fillAcquisitionsData(void* ptr_acq, double v);
	void cSTIR_fillAcquisitionsDataFromAcquisitionsData
		(void* ptr_acq, const void * ptr_from);

	// Reconstruction methods
	void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_runReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_updateReconstruction(void* ptr_r, void* ptr_i);

	// Ojective function methods
	void* cSTIR_setupObjectiveFunction(void* ptr_r, void* ptr_i);
	void*	cSTIR_subsetSensitivity(void* ptr_f, int subset);
	void* cSTIR_objectiveFunctionValue(void* ptr_f, void* ptr_i);
	void* cSTIR_objectiveFunctionGradient
		(void* ptr_f, void* ptr_i, int subset);
	void* cSTIR_objectiveFunctionGradientNotDivided
		(void* ptr_f, void* ptr_i, int subset);

	// Image methods
	void cSTIR_getImageDimensions(const void* ptr, PTR_INT ptr_data);
	void cSTIR_getImageData(const void* ptr, PTR_DOUBLE ptr_data);
	void cSTIR_setImageData(const void* ptr_im, PTR_DOUBLE ptr_data);
	void* cSTIR_voxels3DF(int nx, int ny, int nz,
		double sx, double sy, double sz, double x, double y, double z);
	void* cSTIR_imageFromVoxels(void* ptr_v);
	void* cSTIR_imageFromImage(void* ptr_v);
	void* cSTIR_imageFromAcquisitionData(void* ptr_ad);
	void cSTIR_fillImage(void* ptr_i, double v);
	void* cSTIR_addShape(void* ptr_i, void* ptr_s, float v);
	void* cSTIR_imagesDifference(void* first, void* second, int rimsize);

	// TextWriter methods
	void* newTextPrinter(const char* stream);
	void* newTextWriter(const char* stream);
	void openChannel(int channel, void* ptr_w);
	void closeChannel(int channel, void* ptr_w);
	void setWriter(void* ptr_w, int channel);
	void resetWriter();
	void deleteTextPrinter(void* ptr);
	void deleteTextWriter(void* ptr_w);

#ifndef CSTIR_FOR_MATLAB
}
#endif

#endif
