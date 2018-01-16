/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

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

	// ListmodeToSinogram methods
	void* cSTIR_setListmodeToSinogramsInterval
		(void* ptr_acq, PTR_FLOAT ptr_data);
	void* cSTIR_setupListmodeToSinogramsConverter(void* ptr);
	void* cSTIR_convertListmodeToSinograms(void* ptr);

	// Data processor methods
	void* cSTIR_applyImageDataProcessor(const void* ptr_p, void* ptr_d);

	// Acquisition model methods
	void* cSTIR_createPETAcquisitionSensitivityModel
		(const void* ptr_src, const char* src);
	void* cSTIR_setupAcquisitionSensitivityModel(void* ptr_sm, void* ptr_ad);
	void* cSTIR_applyAcquisitionSensitivityModel(void* ptr_sm, void* ptr_ad);
	void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im);
	void* cSTIR_acquisitionModelFwd(void* ptr_am, void* ptr_im);
	void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad);

	// Acquisition data methods
	void* cSTIR_setAcquisitionsStorageScheme(const char* scheme);
	void* cSTIR_acquisitionsDataFromTemplate(void* ptr_t);
	void* cSTIR_acquisitionsDataFromScannerInfo
		(const char* scanner, int span, int max_ring_diff, int view_mash_factor);
	void* cSTIR_getAcquisitionsDimensions(const void* ptr_acq, PTR_INT ptr_dim);
	void* cSTIR_getAcquisitionsData(const void* ptr_acq, PTR_FLOAT ptr_data);
	void* cSTIR_setAcquisitionsData(void* ptr_acq, PTR_FLOAT ptr_data);
	void* cSTIR_fillAcquisitionsData(void* ptr_acq, float v);
	void* cSTIR_fillAcquisitionsDataFromAcquisitionsData
		(void* ptr_acq, const void * ptr_from);
	void* cSTIR_writeAcquisitionData(void* ptr_acq, const char* filename);

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

	// Prior methods
	void* cSTIR_priorGradient(void* ptr_p, void* ptr_i);

	// Image methods
	void* cSTIR_getImageDimensions(const void* ptr, PTR_INT ptr_data);
	void* cSTIR_getImageData(const void* ptr, PTR_FLOAT ptr_data);
	void* cSTIR_setImageData(const void* ptr_im, PTR_FLOAT ptr_data);
	void* cSTIR_voxels3DF(int nx, int ny, int nz,
		float sx, float sy, float sz, float x, float y, float z);
	void* cSTIR_imageFromVoxels(void* ptr_v);
	void* cSTIR_imageFromImage(void* ptr_v);
	void* cSTIR_imageFromAcquisitionData(void* ptr_ad);
	void* cSTIR_fillImage(void* ptr_i, float v);
	void* cSTIR_addShape(void* ptr_i, void* ptr_s, float v);
	//void* cSTIR_imagesDifference(void* first, void* second, int rimsize);
	void* cSTIR_writeImage(void* ptr_i, const char* filename); 

	// Data container methods
	void* cSTIR_norm(const void* ptr_x);
	void*	cSTIR_dot(const void* ptr_x, const void* ptr_y);
	void* cSTIR_mult(float a, const void* ptr_x);
	void* cSTIR_axpby(float a, const void* ptr_x, float b, const void* ptr_y);

	// TextWriter methods
	void* newTextPrinter(const char* stream);
	void* newTextWriter(const char* stream);
	void openChannel(int channel, void* ptr_w);
	void closeChannel(int channel, void* ptr_w);
	void* deleteTextPrinter(void* ptr);
	void* deleteTextWriter(void* ptr_w);

#ifndef CSTIR_FOR_MATLAB
}
#endif

#endif
