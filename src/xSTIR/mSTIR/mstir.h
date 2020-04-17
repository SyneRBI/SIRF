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
EXPORTED_FUNCTION  void* mSetParameter (void* ptr, const char* obj, const char* name, const void* value);
EXPORTED_FUNCTION 	void* mParameter(const void* ptr, const char* obj, const char* name);
EXPORTED_FUNCTION     void* mSTIR_setVerbosity(const int verbosity_ptr);
EXPORTED_FUNCTION 	void* mSTIR_newObject(const char* name);
EXPORTED_FUNCTION 	void* mSTIR_objectFromFile(const char* name, const char* filename);
EXPORTED_FUNCTION 	void* mSTIR_setParameter (void* ptr, const char* obj, const char* name, const void* value);
EXPORTED_FUNCTION 	void* mSTIR_parameter(const void* ptr, const char* obj, const char* name);
EXPORTED_FUNCTION 	void* mSTIR_setListmodeToSinogramsInterval (void* ptr_acq, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_setListmodeToSinogramsFlag (void* ptr_lm2s, const char* flag, int v);
EXPORTED_FUNCTION 	void* mSTIR_setupListmodeToSinogramsConverter(void* ptr);
EXPORTED_FUNCTION 	void* mSTIR_convertListmodeToSinograms(void* ptr);
EXPORTED_FUNCTION 	void* mSTIR_computeRandoms(void* ptr);
EXPORTED_FUNCTION     void* mSTIR_lm_prompt_rate_exceeds_threshold(void* ptr, const float threshold);
EXPORTED_FUNCTION 	void* mSTIR_setupImageDataProcessor(const void* ptr_p, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_applyImageDataProcessor(const void* ptr_p, void* ptr_d);
EXPORTED_FUNCTION 	void* mSTIR_createPETAcquisitionSensitivityModel (const void* ptr_src, const char* src);
EXPORTED_FUNCTION 	void* mSTIR_createPETAttenuationModel (const void* ptr_img, const void* ptr_am);
EXPORTED_FUNCTION 	void* mSTIR_chainPETAcquisitionSensitivityModels (const void* ptr_first, const void* ptr_second);
EXPORTED_FUNCTION 	void* mSTIR_setupAcquisitionSensitivityModel(void* ptr_sm, void* ptr_ad);
EXPORTED_FUNCTION 	void* mSTIR_applyAcquisitionSensitivityModel (void* ptr_sm, void* ptr_ad, const char* job);
EXPORTED_FUNCTION 	void* mSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelFwd(void* ptr_am, void* ptr_im,  int subset_num, int num_subsets);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelFwdReplace (void* ptr_am, void* ptr_im, int subset_num, int num_subsets, void* ptr_ad);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad, int subset_num, int num_subsets);
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionDataStorageScheme();
EXPORTED_FUNCTION 	void* mSTIR_setAcquisitionDataStorageScheme(const char* scheme);
EXPORTED_FUNCTION 	void* mSTIR_acquisitionDataFromTemplate(void* ptr_t);
EXPORTED_FUNCTION 	void* mSTIR_cloneAcquisitionData(void* ptr_ad);
EXPORTED_FUNCTION 	void* mSTIR_rebinnedAcquisitionData(void* ptr_t, const int num_segments_to_combine, const int num_views_to_combine, const int num_tang_poss_to_trim, const bool do_normalisation, const int max_in_segment_num_to_process );
EXPORTED_FUNCTION 	void* mSTIR_acquisitionDataFromScannerInfo (const char* scanner, int span, int max_ring_diff, int view_mash_factor);
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionDataDimensions(const void* ptr_acq, PTR_INT ptr_dim);
EXPORTED_FUNCTION 	void* mSTIR_getAcquisitionData(const void* ptr_acq, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_setAcquisitionData(void* ptr_acq, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_fillAcquisitionData(void* ptr_acq, float v);
EXPORTED_FUNCTION 	void* mSTIR_fillAcquisitionDataFromAcquisitionData (void* ptr_acq, const void * ptr_from);
EXPORTED_FUNCTION 	void* mSTIR_writeAcquisitionData(void* ptr_acq, const char* filename);
EXPORTED_FUNCTION 	void* mSTIR_setupFBP2DReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_runFBP2DReconstruction(void* ptr_r);
EXPORTED_FUNCTION 	void* mSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_runReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_updateReconstruction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_setupObjectiveFunction(void* ptr_r, void* ptr_i);
EXPORTED_FUNCTION 	void*	mSTIR_subsetSensitivity(void* ptr_f, int subset);
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionValue(void* ptr_f, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionGradient (void* ptr_f, void* ptr_i, int subset);
EXPORTED_FUNCTION 	void* mSTIR_objectiveFunctionGradientNotDivided (void* ptr_f, void* ptr_i, int subset);
EXPORTED_FUNCTION 	void* mSTIR_setupPrior(void* ptr_p, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_priorGradient(void* ptr_p, void* ptr_i);
EXPORTED_FUNCTION 	void* mSTIR_PLSPriorGradient(void* ptr_p, int dir);
EXPORTED_FUNCTION 	void* mSTIR_getImageDimensions(const void* ptr, PTR_INT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_getImageVoxelSizes(const void* ptr_im, PTR_FLOAT ptr_vs);
EXPORTED_FUNCTION 	void* mSTIR_getImageTransformMatrix(const void* ptr_im, PTR_FLOAT ptr_md);
EXPORTED_FUNCTION 	void* mSTIR_getImageData(const void* ptr, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_setImageData(void* ptr_im, PTR_FLOAT ptr_data);
EXPORTED_FUNCTION 	void* mSTIR_setImageDataFromImage(void* ptr_im, const void* ptr_src);
EXPORTED_FUNCTION 	void* mSTIR_voxels3DF(int nx, int ny, int nz, float sx, float sy, float sz, float x, float y, float z);
EXPORTED_FUNCTION 	void* mSTIR_imageFromVoxels(void* ptr_v);
EXPORTED_FUNCTION 	void* mSTIR_imageFromImage(void* ptr_v);
EXPORTED_FUNCTION 	void* mSTIR_imageFromImageData(void* ptr_v);
EXPORTED_FUNCTION 	void* mSTIR_imageFromAcquisitionData(void* ptr_ad);
EXPORTED_FUNCTION 	void* mSTIR_imageFromAcquisitionDataAndNxNy(void* ptr_ad, int nx, int ny);
EXPORTED_FUNCTION 	void* mSTIR_fillImage(void* ptr_i, float v);
EXPORTED_FUNCTION 	void* mSTIR_addShape(void* ptr_i, void* ptr_s, float v);
EXPORTED_FUNCTION 	void* mSTIR_writeImage(void* ptr_i, const char* filename); 
EXPORTED_FUNCTION     void* mSTIR_ImageData_zoom_image(void* ptr_im, const PTR_FLOAT zooms_ptr_raw, const PTR_FLOAT offsets_in_mm_ptr_raw, const PTR_INT new_sizes_ptr_raw, const char * const zoom_options);
EXPORTED_FUNCTION     void* mSTIR_ImageData_move_to_scanner_centre(void* im_ptr, const void* acq_data_ptr);
EXPORTED_FUNCTION 	void* mNewTextPrinter(const char* stream);
EXPORTED_FUNCTION 	void* mNewTextWriter(const char* stream);
EXPORTED_FUNCTION 	void mOpenChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION 	void mCloseChannel(int channel, void* ptr_w);
EXPORTED_FUNCTION 	void* mDeleteTextPrinter(void* ptr);
EXPORTED_FUNCTION 	void* mDeleteTextWriter(void* ptr_w);
#ifndef CSTIR_FOR_MATLAB
}
#endif
EXPORTED_FUNCTION void* mNewMexPrinter();
EXPORTED_FUNCTION void* mDeleteMexPrinter(void* ptr);

#endif
