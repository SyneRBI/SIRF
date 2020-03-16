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

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/STIR/stir_types.h"
#include "sirf/STIR/cstir_p.h"
#include "sirf/STIR/stir_x.h"
#include "stir/ImagingModality.h"
#include "stir/Verbosity.h"

using namespace stir;
using namespace sirf;

#define NEW_OBJECT_HANDLE(T) new ObjectHandle<T >(shared_ptr<T >(new T))
#define SPTR_FROM_HANDLE(Object, X, H) \
	shared_ptr<Object> X; getObjectSptrFromHandle<Object>(H, X);

static void*
unknownObject(const char* obj, const char* name, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "unknown ";
	error += obj;
	error += " '";
	error += name;
	error += "'";
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

template<class Method>
void*
cSTIR_newReconstructionMethod(const char* par_file)
{
	try {
		if (strlen(par_file) > 0) {
			shared_ptr<Reconstruction<Image3DF> > sptr(new Method(par_file));
			return newObjectHandle(sptr);
		}
		else {
			shared_ptr<Reconstruction<Image3DF> > sptr(new Method);
			return newObjectHandle(sptr);
		}
	}
	CATCH;
}

extern "C"
void* cSTIR_setVerbosity(const int verbosity)
{
    stir::Verbosity::set(verbosity);
    return new DataHandle;
}

extern "C"
void* cSTIR_newObject(const char* name)
{
	try {
		if (boost::iequals(name, "FBP2D"))
			return NEW_OBJECT_HANDLE(xSTIR_FBP2DReconstruction);
		if (boost::iequals(name, "ListmodeToSinograms"))
			return NEW_OBJECT_HANDLE(ListmodeToSinograms);
		if (boost::iequals(name,
			"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
			return NEW_OBJECT_HANDLE
			(xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF);
		if (boost::iequals(name, "AcqModUsingMatrix"))
			return NEW_OBJECT_HANDLE(AcqModUsingMatrix3DF);
#ifdef STIR_WITH_NIFTYPET_PROJECTOR
        if (boost::iequals(name, "AcqModUsingNiftyPET"))
            return NEW_OBJECT_HANDLE(AcqModUsingNiftyPET3DF);
#endif
		if (boost::iequals(name, "RayTracingMatrix"))
			return NEW_OBJECT_HANDLE(RayTracingMatrix);
		if (boost::iequals(name, "QuadraticPrior"))
			return NEW_OBJECT_HANDLE(QuadPrior3DF);
		if (boost::iequals(name, "PLSPrior"))
			return NEW_OBJECT_HANDLE(PLSPrior<float>);
		if (boost::iequals(name, "TruncateToCylindricalFOVImageProcessor"))
			return NEW_OBJECT_HANDLE(CylindricFilter3DF);
		if (boost::iequals(name, "EllipsoidalCylinder"))
			return NEW_OBJECT_HANDLE(EllipsoidalCylinder);
		if (boost::iequals(name, "SeparableGaussianImageFilter"))
			return NEW_OBJECT_HANDLE(xSTIR_SeparableGaussianImageFilter);
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void* cSTIR_setParameter
(void* ptr_s, const char* obj, const char* name, const void* ptr_v)
{
	try {
		CAST_PTR(DataHandle, hs, ptr_s);
		CAST_PTR(DataHandle, hv, ptr_v);
		if (boost::iequals(obj, "ListmodeToSinograms"))
			return cSTIR_setListmodeToSinogramsParameter(ptr_s, name, ptr_v);
		else if (boost::iequals(obj, "SeparableGaussianImageFilter"))
			return cSTIR_setSeparableGaussianImageFilterParameter(ptr_s, name, ptr_v);
		else if (boost::iequals(obj, "Shape"))
			return cSTIR_setShapeParameter(ptr_s, name, ptr_v);
		else if (boost::iequals(obj, "EllipsoidalCylinder"))
			return cSTIR_setEllipsoidalCylinderParameter(hs, name, hv);
		else if (boost::iequals(obj, "TruncateToCylindricalFOVImageProcessor"))
			return cSTIR_setTruncateToCylindricalFOVImageProcessorParameter
			(hs, name, hv);
		else if (boost::iequals(obj, "AcquisitionModel"))
			return cSTIR_setAcquisitionModelParameter(hs, name, hv);
		else if (boost::iequals(obj, "AcqModUsingMatrix"))
			return cSTIR_setAcqModUsingMatrixParameter(hs, name, hv);
#ifdef STIR_WITH_NIFTYPET_PROJECTOR
        else if (boost::iequals(obj, "AcqModUsingNiftyPET"))
            return cSTIR_setAcqModUsingNiftyPETParameter(hs, name, hv);
#endif
		else if (boost::iequals(obj, "RayTracingMatrix"))
			return cSTIR_setRayTracingMatrixParameter(hs, name, hv);
		else if (boost::iequals(obj, "GeneralisedPrior"))
			return cSTIR_setGeneralisedPriorParameter(hs, name, hv);
		else if (boost::iequals(obj, "QuadraticPrior"))
			return cSTIR_setQuadraticPriorParameter(hs, name, hv);
		else if (boost::iequals(obj, "PLSPrior"))
			return cSTIR_setPLSPriorParameter(hs, name, hv);
		else if (boost::iequals(obj, "GeneralisedObjectiveFunction"))
			return cSTIR_setGeneralisedObjectiveFunctionParameter(hs, name, hv);
		else if (boost::iequals(obj, "PoissonLogLikelihoodWithLinearModelForMean"))
			return cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanParameter
			(hs, name, hv);
		else if (boost::iequals(obj,
			"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
			return
			cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
			(hs, name, hv);
		else if (boost::iequals(obj, "Reconstruction"))
			return cSTIR_setReconstructionParameter(hs, name, hv);
		else if (boost::iequals(obj, "IterativeReconstruction"))
			return cSTIR_setIterativeReconstructionParameter(hs, name, hv);
		else if (boost::iequals(obj, "OSMAPOSL"))
			return cSTIR_setOSMAPOSLParameter(hs, name, hv);
#ifdef USE_HKEM
		else if (boost::iequals(obj, "KOSMAPOSL"))
			return cSTIR_setKOSMAPOSLParameter(hs, name, hv);
#endif
		else if (boost::iequals(obj, "OSSPS"))
			return cSTIR_setOSSPSParameter(hs, name, hv);
		else if (boost::iequals(obj, "FBP2D"))
			return cSTIR_setFBP2DParameter(hs, name, hv);
		else
			return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void* cSTIR_parameter(const void* ptr, const char* obj, const char* name) 
{
	try {
		CAST_PTR(DataHandle, handle, ptr);
		if (boost::iequals(obj, "Shape"))
			return cSTIR_shapeParameter(handle, name);
		else if (boost::iequals(obj, "EllipsoidalCylinder"))
			return cSTIR_ellipsoidalCylinderParameter(handle, name);
		else if (boost::iequals(obj, "TruncateToCylindricalFOVImageProcessor"))
			return cSTIR_truncateToCylindricalFOVImageProcessorParameter
			(handle, name);
		else if (boost::iequals(obj, "RayTracingMatrix"))
			return cSTIR_rayTracingMatrixParameter(handle, name);
		else if (boost::iequals(obj, "AcqModUsingMatrix"))
			return cSTIR_acqModUsingMatrixParameter(handle, name);
		else if (boost::iequals(obj, "GeneralisedPrior"))
			return cSTIR_generalisedPriorParameter(handle, name);
		else if (boost::iequals(obj, "PLSPrior"))
			return cSTIR_PLSPriorParameter(handle, name);
		else if (boost::iequals(obj, "GeneralisedObjectiveFunction"))
			return cSTIR_generalisedObjectiveFunctionParameter(handle, name);
		else if (boost::iequals(obj,
			"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
			return
			cSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
			(handle, name);
		else if (boost::iequals(obj, "IterativeReconstruction"))
			return cSTIR_iterativeReconstructionParameter(handle, name);
		else if (boost::iequals(obj, "OSMAPOSL"))
			return cSTIR_OSMAPOSLParameter(handle, name);
		else if (boost::iequals(obj, "OSSPS"))
			return cSTIR_OSSPSParameter(handle, name);
		else if (boost::iequals(obj, "FBP2D"))
			return cSTIR_FBP2DParameter(handle, name);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void* cSTIR_objectFromFile(const char* name, const char* filename)
{
	try {
		if (boost::iequals(name, "OSMAPOSLReconstruction"))
			return cSTIR_newReconstructionMethod
			<OSMAPOSLReconstruction<Image3DF> >
			(filename);
#ifdef USE_HKEM
		if (boost::iequals(name, "KOSMAPOSLReconstruction"))
			return cSTIR_newReconstructionMethod
			<KOSMAPOSLReconstruction<Image3DF> >
			(filename);
#endif
		if (boost::iequals(name, "OSSPSReconstruction"))
			return cSTIR_newReconstructionMethod
			<OSSPSReconstruction<Image3DF> >
			(filename);
		if (boost::iequals(name, "Image")) {
			shared_ptr<STIRImageData> sptr(new STIRImageData(filename));
			return newObjectHandle(sptr);
		}
		if (boost::iequals(name, "AcquisitionData")) {
			shared_ptr<PETAcquisitionData> 
				sptr(new PETAcquisitionDataInFile(filename));
			return newObjectHandle(sptr);
		}
		if (boost::iequals(name, "ListmodeToSinograms")) {
			shared_ptr<ListmodeToSinograms>
				sptr(new ListmodeToSinograms(filename));
			return newObjectHandle(sptr);
		}
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

extern "C"
void* cSTIR_setListmodeToSinogramsInterval(void* ptr_lm2s, size_t ptr_data)
{
	try {
		ListmodeToSinograms& lm2s = 
			objectFromHandle<ListmodeToSinograms>(ptr_lm2s);
		float *data = (float *)ptr_data;
		lm2s.set_time_interval((double)data[0], (double)data[1]);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setListmodeToSinogramsFlag(void* ptr_lm2s, const char* flag, int v)
{
	try {
		DataHandle* handle = new DataHandle;
		ListmodeToSinograms& lm2s =
			objectFromHandle<ListmodeToSinograms>(ptr_lm2s);
		int err = lm2s.set_flag(flag, (bool)v);
		if (err) {
			std::string err_msg;
			err_msg = "ListmodeToSinogram does not have this flag: ";
			err_msg += flag;
			ExecutionStatus status(err_msg.c_str(), __FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupListmodeToSinogramsConverter(void* ptr)
{
	try {
		ListmodeToSinograms& lm2s = objectFromHandle<ListmodeToSinograms>(ptr);
		DataHandle* handle = new DataHandle;
		if (lm2s.set_up()) {
			ExecutionStatus status
				("cSTIR_setupListmodeToSinogramConverter failed", 
					__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_convertListmodeToSinograms(void* ptr)
{
	try {
		ListmodeToSinograms& lm2s = objectFromHandle<ListmodeToSinograms>(ptr);
		lm2s.process_data();
		return newObjectHandle(lm2s.get_output());
	}
	CATCH;
}

extern "C"
void* cSTIR_computeRandoms(void* ptr)
{
	try {
		ListmodeToSinograms& lm2s = objectFromHandle<ListmodeToSinograms>(ptr);
		if (lm2s.estimate_randoms()) {
			ExecutionStatus status
				("cSTIR_computeRandoms failed", __FILE__, __LINE__);
			DataHandle* handle = new DataHandle;
			handle->set(0, &status);
			return handle;
		}
		return newObjectHandle(lm2s.get_randoms_sptr());
	}
	CATCH;
}

extern "C"
void* cSTIR_lm_prompt_rate_exceeds_threshold(const void * ptr, const float threshold)
{
    try {
        ListmodeToSinograms& lm2s = objectFromHandle<ListmodeToSinograms>(ptr);
        return dataHandle<float>(lm2s.get_time_at_which_prompt_rate_exceeds_threshold(threshold));
    }
    CATCH
}
extern "C"
void* cSTIR_setupImageDataProcessor(const void* ptr_p, void* ptr_i)
{
	try {
		DataProcessor<Image3DF>& processor =
			objectFromHandle<DataProcessor<Image3DF> >(ptr_p);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		processor.set_up(image);
		return (void*) new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_applyImageDataProcessor(const void* ptr_p, void* ptr_i)
{
	try {
		DataProcessor<Image3DF>& processor =
			objectFromHandle<DataProcessor<Image3DF> >(ptr_p);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		processor.apply(image);
		return (void*) new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_createPETAcquisitionSensitivityModel
	(const void* ptr_src, const char* src)
{
	try {
		shared_ptr<PETAcquisitionSensitivityModel> sptr;
		if (boost::iequals(src, "s")) {
			PETAcquisitionData& ad = objectFromHandle<PETAcquisitionData>(ptr_src);
			sptr.reset(new PETAcquisitionSensitivityModel(ad));
		}
		else if (boost::iequals(src, "n")) {
			CAST_PTR(DataHandle, h, ptr_src);
			sptr.reset(new PETAcquisitionSensitivityModel(charDataFromDataHandle(h)));
		}
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_createPETAttenuationModel(const void* ptr_img, const void* ptr_am)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_img);
		PETAcquisitionModel& am = objectFromHandle<PETAcquisitionModel>(ptr_am);
		shared_ptr<PETAcquisitionSensitivityModel> 
			sptr(new PETAttenuationModel(id, am));
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_chainPETAcquisitionSensitivityModels
(const void* ptr_first, const void* ptr_second)
{
	try {
		PETAcquisitionSensitivityModel& first =
			objectFromHandle<PETAcquisitionSensitivityModel>(ptr_first);
		PETAcquisitionSensitivityModel& second =
			objectFromHandle<PETAcquisitionSensitivityModel>(ptr_second);
		shared_ptr<PETAcquisitionSensitivityModel> 
			sptr(new PETAcquisitionSensitivityModel(first, second));
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_setupAcquisitionSensitivityModel(void* ptr_sm, void* ptr_ad)
{
	try {
		PETAcquisitionSensitivityModel& sm = 
			objectFromHandle<PETAcquisitionSensitivityModel>(ptr_sm);
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_ad);
		Succeeded s = sm.set_up(sptr_ad->data()->get_proj_data_info_sptr());
		DataHandle* handle = new DataHandle;
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_acquisitionModelSetup failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_applyAcquisitionSensitivityModel
(void* ptr_sm, void* ptr_ad, const char* job)
{
	try {
		void* handle = new DataHandle;
		PETAcquisitionSensitivityModel& sm =
			objectFromHandle<PETAcquisitionSensitivityModel>(ptr_sm);
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_ad);
		if (boost::iequals(job, "unnormalise"))
			sm.unnormalise(*sptr_ad);
		else if (boost::iequals(job, "normalise"))
			sm.normalise(*sptr_ad);
		else if (boost::iequals(job, "fwd"))
			handle = newObjectHandle(sm.forward(*sptr_ad));
		else if (boost::iequals(job, "inv"))
			handle = newObjectHandle(sm.invert(*sptr_ad));
		return handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im)
{
	try {
		//writeText("setting up acquisition model\n");
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_dt, ptr_dt);
		SPTR_FROM_HANDLE(STIRImageData, sptr_id, ptr_im);
		//std::cout << "setting up acquisition model...\n";
		Succeeded s = am.set_up(sptr_dt, sptr_id);
		DataHandle* handle = new DataHandle;
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_acquisitionModelSetup failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionModelFwd
(void* ptr_am, void* ptr_im, int subset_num, int num_subsets)
{
	try {
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		return newObjectHandle(am.forward(id, subset_num, num_subsets));
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionModelFwdReplace
(void* ptr_am, void* ptr_im, int subset_num, int num_subsets, void* ptr_ad)
{
	try {
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		PETAcquisitionData& ad = objectFromHandle<PETAcquisitionData>(ptr_ad);
		am.forward(ad, id, subset_num, num_subsets);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad, 
	int subset_num, int num_subsets)
{
	try {
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		PETAcquisitionData& ad = objectFromHandle<PETAcquisitionData>(ptr_ad);
		return newObjectHandle(am.backward(ad, subset_num, num_subsets));
	}
	CATCH;
}

extern "C"
void*
cSTIR_setAcquisitionDataStorageScheme(const char* scheme)
{ 
	try {
		if (scheme[0] == 'f' || strcmp(scheme, "default") == 0)
			PETAcquisitionDataInFile::set_as_template();
		else
			PETAcquisitionDataInMemory::set_as_template();
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void*
cSTIR_getAcquisitionDataStorageScheme()
{
	return charDataHandleFromCharData
		(PETAcquisitionData::storage_scheme().c_str());
}

extern "C"
void* cSTIR_acquisitionDataFromTemplate(void* ptr_t)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_t, ptr_t);
		shared_ptr<PETAcquisitionData> sptr(sptr_t->new_acquisition_data());
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_cloneAcquisitionData(void* ptr_ad)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_ad);
		shared_ptr<PETAcquisitionData> sptr(sptr_ad->clone());
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_rebinnedAcquisitionData(void* ptr_t, 
const int num_segments_to_combine,
const int num_views_to_combine,
const int num_tang_poss_to_trim,
const bool do_normalisation,
const int max_in_segment_num_to_process
)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_t, ptr_t);
		shared_ptr<PETAcquisitionData> sptr =
			sptr_t->single_slice_rebinned_data(
			num_segments_to_combine,
			num_views_to_combine,
			num_tang_poss_to_trim,
			do_normalisation,
			max_in_segment_num_to_process
			);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionDataFromScannerInfo
(const char* scanner, int span, int max_ring_diff, int view_mash_factor)
{
	try{
		shared_ptr<ExamInfo> sptr_ei(new ExamInfo());
        sptr_ei->imaging_modality = ImagingModality::PT;
		stir::shared_ptr<stir::ProjDataInfo> sptr_pdi =
			PETAcquisitionData::proj_data_info_from_scanner
			(scanner, span, max_ring_diff, view_mash_factor);
		PETAcquisitionDataInFile::init();
		stir::shared_ptr<PETAcquisitionData> sptr_t =
			PETAcquisitionData::storage_template();
		stir::shared_ptr<PETAcquisitionData> sptr(sptr_t->same_acquisition_data
			(sptr_ei, sptr_pdi));
		sptr->fill(0.0f);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_getAcquisitionDataDimensions(const void* ptr_acq, size_t ptr_dim)
{
	try {
		int* dim = (int*)ptr_dim;
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		dim[0] = sptr_ad->get_num_tangential_poss();
		dim[1] = sptr_ad->get_num_views();
		dim[2] = sptr_ad->get_num_sinograms();
		dim[3] = sptr_ad->get_num_TOF_bins();
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getAcquisitionData(const void* ptr_acq, size_t ptr_data)
{
	try {
		float* data = (float*)ptr_data;
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		sptr_ad->copy_to(data);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_fillAcquisitionData(void* ptr_acq, float v)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		sptr_ad->fill(v);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_fillAcquisitionDataFromAcquisitionData
(void* ptr_acq, const void* ptr_from)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_from, ptr_from);
		sptr_ad->fill(*sptr_from);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setAcquisitionData(void* ptr_acq, size_t ptr_data)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		float *data = (float *)ptr_data;
		sptr_ad->fill_from(data);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_writeAcquisitionData(void* ptr_acq, const char* filename)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		sptr_ad->write(filename);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupFBP2DReconstruction(void* ptr_r, void* ptr_i)
{
	try {
		DataHandle* handle = new DataHandle;
		xSTIR_FBP2DReconstruction& recon =
			objectFromHandle< xSTIR_FBP2DReconstruction >(ptr_r);
		SPTR_FROM_HANDLE(STIRImageData, sptr_id, ptr_i);
		if (recon.set_up(sptr_id) != Succeeded::yes) {
			ExecutionStatus status("cSTIR_setupFBP2DReconstruction failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_runFBP2DReconstruction(void* ptr_r)
{
	try {
		DataHandle* handle = new DataHandle;
		xSTIR_FBP2DReconstruction& recon =
			objectFromHandle< xSTIR_FBP2DReconstruction >(ptr_r);
		if (recon.process() != Succeeded::yes) {
			ExecutionStatus status("cSTIR_FBP2DReconstruction failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i)
{
	try {
		//std::cout << "in cSTIR_setupReconstruction...\n";
		DataHandle* handle = new DataHandle;
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		sptrImage3DF sptr_image = id.data_sptr();
		xSTIR_IterativeReconstruction3DF& recon =
			objectFromHandle<xSTIR_IterativeReconstruction3DF>(ptr_r);
		Succeeded s = Succeeded::no;
		if (!recon.post_process()) {
			s = recon.setup(sptr_image);
			recon.subiteration() = recon.get_start_subiteration_num();
		}
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_setupReconstruction failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_runReconstruction(void* ptr_r, void* ptr_i) 
{
	try {
		DataHandle* handle = new DataHandle;
		Reconstruction<Image3DF>& recon =
			objectFromHandle< Reconstruction<Image3DF> >(ptr_r);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		sptrImage3DF sptr_image = id.data_sptr();
		if (recon.reconstruct(sptr_image) != Succeeded::yes) {
			ExecutionStatus status("cSTIR_reconstruct failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void* cSTIR_updateReconstruction(void* ptr_r, void* ptr_i)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		xSTIR_IterativeReconstruction3DF& recon =
			objectFromHandle<xSTIR_IterativeReconstruction3DF>(ptr_r);
		recon.update(image);
		return (void*) new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupObjectiveFunction(void* ptr_r, void* ptr_i)
{
	try {
		DataHandle* handle = new DataHandle;
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		sptrImage3DF sptr_image = id.data_sptr();
		xSTIR_GeneralisedObjectiveFunction3DF& obj_fun =
			objectFromHandle<xSTIR_GeneralisedObjectiveFunction3DF>(ptr_r);
		Succeeded s = Succeeded::no;
		if (!obj_fun.post_process())
			s = obj_fun.set_up(sptr_image);
		if (s != Succeeded::yes) {
			ExecutionStatus status("cSTIR_setupObjectiveFunction failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cSTIR_objectiveFunctionValue(void* ptr_f, void* ptr_i)
{
	try {
		ObjectiveFunction3DF& fun = objectFromHandle< ObjectiveFunction3DF>(ptr_f);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		float v = (float)fun.compute_objective_function(image);
		return dataHandle<float>(v);
	}
	CATCH;
}

extern "C"
void*
cSTIR_subsetSensitivity(void* ptr_f, int subset)
{
	try {
		PoissonLogLhLinModMean3DF& fun =
			objectFromHandle<PoissonLogLhLinModMean3DF>(ptr_f);
		const Image3DF& s = fun.get_subset_sensitivity(subset);
		STIRImageData* ptr_id = new STIRImageData(s);
		shared_ptr<STIRImageData> sptr(ptr_id);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_objectiveFunctionGradient(void* ptr_f, void* ptr_i, int subset)
{
	try {
		ObjectiveFunction3DF& fun = objectFromHandle< ObjectiveFunction3DF>(ptr_f);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		STIRImageData* ptr_id = new STIRImageData(image);
		shared_ptr<STIRImageData> sptr(ptr_id);
		Image3DF& grad = sptr->data();
		if (subset >= 0)
			fun.compute_sub_gradient(grad, image, subset);
		else {
			int nsub = fun.get_num_subsets();
			grad.fill(0.0);
			STIRImageData* ptr_id = new STIRImageData(image);
			shared_ptr<STIRImageData> sptr_sub(ptr_id);
			Image3DF& subgrad = sptr_sub->data();
			for (int sub = 0; sub < nsub; sub++) {
				fun.compute_sub_gradient(subgrad, image, sub);
				grad += subgrad;
			}
		}
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_objectiveFunctionGradientNotDivided(void* ptr_f, void* ptr_i, int subset)
{
	try {
		PoissonLogLhLinModMean3DF& fun = 
			objectFromHandle<PoissonLogLhLinModMean3DF>(ptr_f);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<STIRImageData> sptr(new STIRImageData(image));
		Image3DF& grad = sptr->data();
		fun.compute_sub_gradient_without_penalty_plus_sensitivity
			(grad, image, subset);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_setupPrior(void* ptr_p, void* ptr_i)
{
	try {
		DataHandle* handle = new DataHandle;
		xSTIR_GeneralisedPrior3DF& prior =
			objectFromHandle<xSTIR_GeneralisedPrior3DF>(ptr_p);
		STIRImageData& image = objectFromHandle<STIRImageData>(ptr_i);
		sptrImage3DF sptr_img = image.data_sptr();
		// empty image is a temporary measure for compatibility with old scripts
		// (valid for as long as the argument of prior.set_up() is not used)
		//sptrImage3DF sptr_img(new Voxels3DF);
		prior.set_up(sptr_img);
		if (prior.post_process()){
			ExecutionStatus status("cSTIR_setupPrior failed",
				__FILE__, __LINE__);
			handle->set(0, &status);
		}
		return handle;
	}
	CATCH;
}

extern "C"
void*
cSTIR_priorGradient(void* ptr_p, void* ptr_i)
{
	try {
		Prior3DF& prior = objectFromHandle<Prior3DF>(ptr_p);
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<STIRImageData> sptr(new STIRImageData(image));
		Image3DF& grad = sptr->data();
		prior.compute_gradient(grad, image);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_PLSPriorGradient(void* ptr_p, int dir)
{
	try {
		PLSPrior<float>& prior = objectFromHandle<PLSPrior<float> >(ptr_p);
		sptrImage3DF sptr_im = prior.get_anatomical_grad_sptr(dir);
		shared_ptr<STIRImageData> sptr_id(new STIRImageData(sptr_im));
		return newObjectHandle(sptr_id);
	}
	CATCH;
}

extern "C"
void* cSTIR_voxels3DF
(int nx, int ny, int nz,
float sx, float sy, float sz,
float x, float y, float z)
{
	try {
		shared_ptr<Voxels3DF> sptr(new Voxels3DF(IndexRange3D(0, nz - 1,
			-(ny / 2), -(ny / 2) + ny - 1, -(nx / 2), -(nx / 2) + nx - 1),
			Coord3DF(z, y, x),
			Coord3DF(sz, sy, sx)));
		sptr->fill(0.0);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromVoxels(void* ptr_v)
{
	try {
		Voxels3DF& voxels = objectFromHandle<Voxels3DF>(ptr_v);
		shared_ptr<STIRImageData> sptr(new STIRImageData(voxels));
		return (void*)newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromImageData(void* ptr_v)
{
	try {
		ImageData& id = objectFromHandle<ImageData>(ptr_v);
		shared_ptr<STIRImageData> sptr(new STIRImageData(id));
		return (void*)newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromImage(void* ptr_i)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<STIRImageData> sptr(new STIRImageData(image));
		return (void*)newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_writeImage(void* ptr_i, const char* filename)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		id.write(filename);
		return (void*) new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_ImageData_zoom_image(void* ptr_im, const size_t zooms_ptr_raw, const size_t offsets_in_mm_ptr_raw,
                                 const size_t new_sizes_ptr_raw, const char *const zoom_options)
{
    try {
        STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);

        const float* zooms_ptr         = (const float*)zooms_ptr_raw;
        const float* offsets_in_mm_ptr = (const float*)offsets_in_mm_ptr_raw;
        const  int*  new_sizes_ptr     = (const  int* )new_sizes_ptr_raw;

        Coord3DF zooms(zooms_ptr[0],zooms_ptr[1],zooms_ptr[2]);
        Coord3DF offsets_in_mm(offsets_in_mm_ptr[0],offsets_in_mm_ptr[1],offsets_in_mm_ptr[2]);
        Coord3DI new_sizes(new_sizes_ptr[0],new_sizes_ptr[1],new_sizes_ptr[2]);

        id.zoom_image(zooms, offsets_in_mm, new_sizes, zoom_options);

		return static_cast<void*>(new DataHandle);
	}
	CATCH;
}

extern "C"
void* cSTIR_ImageData_move_to_scanner_centre(void* im_ptr, const void* acq_data_ptr)
{
    try {
        STIRImageData& im = objectFromHandle<STIRImageData>(im_ptr);
        PETAcquisitionData& ad = objectFromHandle<PETAcquisitionData>(acq_data_ptr);
        im.move_to_scanner_centre(ad);

        return static_cast<void*>(new DataHandle);
	}
	CATCH;

}

extern "C"
void* cSTIR_imageFromAcquisitionData(void* ptr_ad)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_ad);
		shared_ptr<STIRImageData> sptr(new STIRImageData(*sptr_ad));
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromAcquisitionDataAndNxNy(void* ptr_ad, int nx, int ny)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_ad);
		STIRImageData id(*sptr_ad);
		int dim[3];
		float vs[3];
		float is[3];
		id.get_dimensions(dim);
		id.get_voxel_sizes(vs);
		for (int i = 0; i < 3; i++)
			is[i] = dim[i] * vs[i];
		int nz = dim[0];
		float vx = is[2] / nx;
		float vy = is[1] / ny;
		float vz = vs[0];
		shared_ptr<Voxels3DF> sptr_v(new Voxels3DF(IndexRange3D(0, nz - 1,
			-(ny / 2), -(ny / 2) + ny - 1, -(nx / 2), -(nx / 2) + nx - 1),
			Coord3DF(0, 0, 0),
			Coord3DF(vz, vy, vx)));
		shared_ptr<STIRImageData> sptr(new STIRImageData(*sptr_v));
		sptr->fill(0.0);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_addShape(void* ptr_i, void* ptr_s, float v)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		sptrVoxels3DF sptr_v((Voxels3DF*)image.clone());
		Voxels3DF& voxels = *sptr_v;
		Shape3D& shape = objectFromHandle<Shape3D>(ptr_s);
		CartesianCoordinate3D<int> num_samples(1, 1, 1);
		voxels.fill(0);
		shape.construct_volume(voxels, num_samples);
		voxels *= v;
		image += voxels;
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_fillImage(void* ptr_i, float v)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_i);
		Image3DF& image = id.data();
		image.fill(v);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getImageDimensions(const void* ptr_im, size_t ptr_dim) 
{
	try {
		int* dim = (int*)ptr_dim;
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		if (id.get_dimensions(dim)) {
				ExecutionStatus status("not a regular image", __FILE__, __LINE__);
				DataHandle* handle = new DataHandle;
				handle->set(0, &status);
				return (void*)handle;
		}
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getImageVoxelSizes(const void* ptr_im, size_t ptr_vs)
{
	try {
		float* vs = (float*)ptr_vs;
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		id.get_voxel_sizes(vs);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getImageTransformMatrix(const void* ptr_im, size_t ptr_md)
{
	try {
        STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		float* data = (float*)ptr_md;
		TransformMatrix3D mx = id.get_geom_info_sptr()->calculate_index_to_physical_point_matrix();
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				data[i + 4 * j] = mx[j][i];
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getImageData(const void* ptr_im, size_t ptr_data)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		float* data = (float*)ptr_data;
		id.get_data(data);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setImageData(void* ptr_im, size_t ptr_data)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		float* data = (float*)ptr_data;
		id.set_data(data);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setImageDataFromImage(void* ptr_im, const void* ptr_src)
{
	try {
		STIRImageData& id = objectFromHandle<STIRImageData>(ptr_im);
		STIRImageData& id_src = objectFromHandle<STIRImageData>(ptr_src);
		Image3DF& data = id.data();
		data = id_src.data();
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* setParameter
(void* ptr_s, const char* obj, const char* name, const void* ptr_v)
{
	return cSTIR_setParameter(ptr_s, obj, name, ptr_v);
}

extern "C"
void* parameter(const void* ptr, const char* obj, const char* name)
{
	return cSTIR_parameter(ptr, obj, name);
}
