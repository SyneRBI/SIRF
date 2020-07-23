/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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

	// Unified parameter exchange methods
	void* setParameter
		(void* ptr, const char* obj, const char* name, const void* value);
	void* parameter(const void* ptr, const char* obj, const char* name);

    // Global
    void* cSTIR_setVerbosity(const int verbosity_ptr);
    void* cSTIR_getVerbosity();
    void* cSTIR_setOMPThreads(const int threads);
    void* cSTIR_getOMPThreads();
	void* cSTIR_useDefaultOMPThreads();
	void* cSTIR_getDefaultOMPThreads();

	// Common STIR Object methods
	void* cSTIR_newObject(const char* name);
	void* cSTIR_objectFromFile(const char* name, const char* filename);
	void* cSTIR_setParameter
		(void* ptr, const char* obj, const char* name, const void* value);
	void* cSTIR_parameter(const void* ptr, const char* obj, const char* name);

	// ListmodeToSinogram methods
	void* cSTIR_setListmodeToSinogramsInterval
		(void* ptr_acq, PTR_FLOAT ptr_data);
	void* cSTIR_setListmodeToSinogramsFlag
		(void* ptr_lm2s, const char* flag, int v);
	void* cSTIR_setupListmodeToSinogramsConverter(void* ptr);
	void* cSTIR_convertListmodeToSinograms(void* ptr);
	void* cSTIR_computeRandoms(void* ptr);
    void* cSTIR_lm_num_prompts_exceeds_threshold(void* ptr, const float threshold);

	// Data processor methods
	void* cSTIR_setupImageDataProcessor(const void* ptr_p, void* ptr_i);
	void* cSTIR_applyImageDataProcessor(const void* ptr_p, void* ptr_d);

	// Acquisition model methods
	void* cSTIR_createPETAcquisitionSensitivityModel
		(const void* ptr_src, const char* src);
	void* cSTIR_createPETAttenuationModel
		(const void* ptr_img, const void* ptr_am);
	void* cSTIR_chainPETAcquisitionSensitivityModels
		(const void* ptr_first, const void* ptr_second);
	void* cSTIR_setupAcquisitionSensitivityModel(void* ptr_sm, void* ptr_ad);
	void* cSTIR_applyAcquisitionSensitivityModel
		(void* ptr_sm, void* ptr_ad, const char* job);
	void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im);
	void* cSTIR_acquisitionModelFwd(void* ptr_am, void* ptr_im, 
		int subset_num, int num_subsets);
	void* cSTIR_acquisitionModelFwdReplace
		(void* ptr_am, void* ptr_im, int subset_num, int num_subsets, void* ptr_ad);
	void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad,
		int subset_num, int num_subsets);

	// Acquisition data methods
	void* cSTIR_getAcquisitionDataStorageScheme();
	void* cSTIR_setAcquisitionDataStorageScheme(const char* scheme);
	void* cSTIR_acquisitionDataFromTemplate(void* ptr_t);
	void* cSTIR_cloneAcquisitionData(void* ptr_ad);
	void* cSTIR_rebinnedAcquisitionData(void* ptr_t,
		const int num_segments_to_combine,
		const int num_views_to_combine,
		const int num_tang_poss_to_trim,
		const bool do_normalisation,
		const int max_in_segment_num_to_process
		);
	void* cSTIR_acquisitionDataFromScannerInfo
		(const char* scanner, int span, int max_ring_diff, int view_mash_factor);
	void* cSTIR_getAcquisitionDataDimensions(const void* ptr_acq, PTR_INT ptr_dim);
	void* cSTIR_getAcquisitionData(const void* ptr_acq, PTR_FLOAT ptr_data);
	void* cSTIR_setAcquisitionData(void* ptr_acq, PTR_FLOAT ptr_data);
	void* cSTIR_fillAcquisitionData(void* ptr_acq, float v);
	void* cSTIR_fillAcquisitionDataFromAcquisitionData
		(void* ptr_acq, const void * ptr_from);
	void* cSTIR_writeAcquisitionData(void* ptr_acq, const char* filename);
	void* cSTIR_get_ProjDataInfo(void* ptr_acq);

	// Reconstruction methods
	void* cSTIR_setupFBP2DReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_runFBP2DReconstruction(void* ptr_r);
	void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_runReconstruction(void* ptr_r, void* ptr_i);
	void* cSTIR_updateReconstruction(void* ptr_r, void* ptr_i);

	// Objective function methods
	void* cSTIR_setupObjectiveFunction(void* ptr_r, void* ptr_i);
	void*	cSTIR_subsetSensitivity(void* ptr_f, int subset);
	void* cSTIR_objectiveFunctionValue(void* ptr_f, void* ptr_i);
	void* cSTIR_objectiveFunctionGradient
		(void* ptr_f, void* ptr_i, int subset);
	void* cSTIR_objectiveFunctionGradientNotDivided
		(void* ptr_f, void* ptr_i, int subset);

	// Prior methods
	void* cSTIR_setupPrior(void* ptr_p, void* ptr_i);
	void* cSTIR_priorGradient(void* ptr_p, void* ptr_i);
	void* cSTIR_PLSPriorGradient(void* ptr_p, int dir);

	// Image methods
	void* cSTIR_getImageDimensions(const void* ptr, PTR_INT ptr_data);
	void* cSTIR_getImageVoxelSizes(const void* ptr_im, PTR_FLOAT ptr_vs);
	void* cSTIR_getImageTransformMatrix(const void* ptr_im, PTR_FLOAT ptr_md);
	void* cSTIR_getImageData(const void* ptr, PTR_FLOAT ptr_data);
	void* cSTIR_setImageData(void* ptr_im, PTR_FLOAT ptr_data);
	void* cSTIR_setImageDataFromImage(void* ptr_im, const void* ptr_src);
	void* cSTIR_voxels3DF(int nx, int ny, int nz,
		float sx, float sy, float sz, float x, float y, float z);
	void* cSTIR_imageFromVoxels(void* ptr_v);
	void* cSTIR_imageFromImage(void* ptr_v);
	void* cSTIR_imageFromImageData(void* ptr_v);
	void* cSTIR_imageFromAcquisitionData(void* ptr_ad);
	void* cSTIR_imageFromAcquisitionDataAndNxNy(void* ptr_ad, int nx, int ny);
	void* cSTIR_fillImage(void* ptr_i, float v);
	void* cSTIR_addShape(void* ptr_i, void* ptr_s, float v, int num_samples_in_each_direction);
	void* cSTIR_writeImage(void* ptr_i, const char* filename); 
    void* cSTIR_writeImage_par(void* ptr_i, const char* filename, const char* par);
    void* cSTIR_ImageData_zoom_image(void* ptr_im,
                                     const PTR_FLOAT zooms_ptr_raw,
                                     const PTR_FLOAT offsets_in_mm_ptr_raw,
                                     const PTR_INT new_sizes_ptr_raw,
                                     const char * const zoom_options);
    void* cSTIR_ImageData_move_to_scanner_centre(void* im_ptr, const void* acq_data_ptr);

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
