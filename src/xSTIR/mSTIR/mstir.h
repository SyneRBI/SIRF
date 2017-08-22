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
#ifndef CSTIR_TO_MATLAB_INTERFACE
#define CSTIR_TO_MATLAB_INTERFACE

#define CSTIR_FOR_MATLAB
#ifdef _WIN32
#define EXPORTED_FUNCTION __declspec(dllexport)
#else
#define EXPORTED_FUNCTION
#endif

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
EXPORTED_FUNCTION 	void* mSTIR_applyDataProcessor(const void* ptr_p, void* ptr_d);
EXPORTED_FUNCTION 	void* mSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelFwd (void* ptr_am, void* ptr_im, const char* datafile);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad);
EXPORTED_FUNCTION 	void* mSTIR_setAcquisitionsStorageScheme(const char* scheme);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionsDataFromTemplate(void* ptr_t);
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionsDimensions(const void* ptr_acq, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionsData(const void* ptr_acq, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_setAcquisitionsData(void* ptr_acq, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_fillAcquisitionsData(void* ptr_acq, float v);
EXPORTED_FUNCTION 	void* mSTIR_fillAcquisitionsDataFromAcquisitionsData (void* ptr_acq, const void * ptr_from);
EXPORTED_FUNCTION 	void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_setupObjectiveFunction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void*	mSTIR_subsetSensitivity(void* ptr_f, int subset);
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionValue(void* ptr_f, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionGradient (void* ptr_f, void* ptr_i, int subset);
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionGradientNotDivided (void* ptr_f, void* ptr_i, int subset);
EXPORTED_FUNCTION 	void* mSTIR_priorGradient(void* ptr_p, void* ptr_i);
EXPORTED_FUNCTION 	void mSTIR_getImageDimensions(const void* ptr, PTR_INT ptr_data);
EXPORTED_FUNCTION 	void mSTIR_getImageData(const void* ptr, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void mSTIR_setImageData(const void* ptr_im, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_voxels3DF(int nx, int ny, int nz, float sx, float sy, float sz, float x, float y, float z);
EXPORTED_FUNCTION 	void* mSTIR_imageFromVoxels(void* ptr_v);
EXPORTED_FUNCTION 	void* mSTIR_imageFromImage(void* ptr_v);
EXPORTED_FUNCTION 	void* mSTIR_imageFromAcquisitionData(void* ptr_ad);
EXPORTED_FUNCTION 	void mSTIR_fillImage(void* ptr_i, float v);
EXPORTED_FUNCTION 	void* mSTIR_addShape(void* ptr_i, void* ptr_s, float v);
//EXPORTED_FUNCTION 	void* mSTIR_imagesDifference(void* first, void* second, int rimsize);
EXPORTED_FUNCTION 	void* mSTIR_writeImage(void* ptr_i, const char* filename); 
EXPORTED_FUNCTION 	void* mNewTextPrinter(const char* stream);
EXPORTED_FUNCTION 	void* mNewTextWriter(const char* stream);
EXPORTED_FUNCTION 	void mOpenChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION 	void mCloseChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION 	void mDeleteTextPrinter(void* ptr);
EXPORTED_FUNCTION 	void mDeleteTextWriter(void* ptr_w);
#ifndef CSTIR_FOR_MATLAB
}
#endif
EXPORTED_FUNCTION void* mNewMexPrinter();
EXPORTED_FUNCTION void mDeleteMexPrinter(void* ptr);

#endif
