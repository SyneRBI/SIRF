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

#include "cstir_shared_ptr.h"
#include "data_handle.h"
#include "cstir_p.h"
#include "stir_types.h"
#include "stir_x.h"
#include "stir_data_containers.h"

using stir::shared_ptr;

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
void* cSTIR_newObject(const char* name)
{
	try {
		if (boost::iequals(name, "ListmodeToSinograms"))
			return newObjectHandle<ListmodeToSinograms>();
		if (boost::iequals(name,
			"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
			return newObjectHandle
			<xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF>();
		if (boost::iequals(name, "AcqModUsingMatrix"))
			return newObjectHandle<AcqModUsingMatrix3DF>();
		if (boost::iequals(name, "RayTracingMatrix"))
			return newObjectHandle<RayTracingMatrix>();
		if (boost::iequals(name, "QuadraticPrior"))
			return newObjectHandle<QuadPrior3DF>();
		if (boost::iequals(name, "TruncateToCylindricalFOVImageProcessor"))
			return newObjectHandle<CylindricFilter3DF>();
		if (boost::iequals(name, "EllipsoidalCylinder"))
			return newObjectHandle<EllipsoidalCylinder>();
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
		else if (boost::iequals(obj, "RayTracingMatrix"))
			return cSTIR_setRayTracingMatrixParameter(hs, name, hv);
		else if (boost::iequals(obj, "GeneralisedPrior"))
			return cSTIR_setGeneralisedPriorParameter(hs, name, hv);
		else if (boost::iequals(obj, "QuadraticPrior"))
			return cSTIR_setQuadraticPriorParameter(hs, name, hv);
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
		else if (boost::iequals(obj, "OSSPS"))
			return cSTIR_setOSSPSParameter(hs, name, hv);
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
		if (boost::iequals(obj, "EllipsoidalCylinder"))
			return cSTIR_ellipsoidalCylinderParameter(handle, name);
		else if (boost::iequals(obj, "TruncateToCylindricalFOVImageProcessor"))
			return cSTIR_truncateToCylindricalFOVImageProcessorParameter
			(handle, name);
		if (boost::iequals(obj, "RayTracingMatrix"))
			return cSTIR_rayTracingMatrixParameter(handle, name);
		else if (boost::iequals(obj, "AcqModUsingMatrix"))
			return cSTIR_acqModUsingMatrixParameter(handle, name);
		if (boost::iequals(obj, "GeneralisedPrior"))
			return cSTIR_generalisedPriorParameter(handle, name);
		if (boost::iequals(obj, "GeneralisedObjectiveFunction"))
			return cSTIR_generalisedObjectiveFunctionParameter(handle, name);
		if (boost::iequals(obj,
			"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
			return
			cSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
			(handle, name);
		if (boost::iequals(obj, "IterativeReconstruction"))
			return cSTIR_iterativeReconstructionParameter(handle, name);
		if (boost::iequals(obj, "OSMAPOSL"))
			return cSTIR_OSMAPOSLParameter(handle, name);
		if (boost::iequals(obj, "OSSPS"))
			return cSTIR_OSSPSParameter(handle, name);
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
		if (boost::iequals(name, "OSSPSReconstruction"))
			return cSTIR_newReconstructionMethod
			<OSSPSReconstruction<Image3DF> >
			(filename);
		if (boost::iequals(name, "Image")) {
			shared_ptr<PETImageData> sptr(new PETImageData(filename));
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
		//return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_computeRandoms(void* ptr)
{
	try {
		ListmodeToSinograms& lm2s = objectFromHandle<ListmodeToSinograms>(ptr);
		lm2s.compute_fan_sums();
		lm2s.compute_singles();
		lm2s.compute_randoms();
		return newObjectHandle(lm2s.get_randoms_sptr());
	}
	CATCH;
}

extern "C"
void* cSTIR_applyImageDataProcessor(const void* ptr_p, void* ptr_i)
{
	try {
		DataProcessor<Image3DF>& processor =
			objectFromHandle<DataProcessor<Image3DF> >(ptr_p);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		else if (boost::iequals(src, "i")) {
			PETImageData& id = objectFromHandle<PETImageData>(ptr_src);
			sptr.reset(new PETAcquisitionSensitivityModel(id));
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
		SPTR_FROM_HANDLE(PETImageData, sptr_id, ptr_im);
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
(void* ptr_am, void* ptr_im)
{
	try {
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		return newObjectHandle(am.forward(id));
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad)
{
	try {
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		PETAcquisitionData& ad = objectFromHandle<PETAcquisitionData>(ptr_ad);
		return newObjectHandle(am.backward(ad));
	}
	CATCH;
}

extern "C"
void*
cSTIR_setAcquisitionsStorageScheme(const char* scheme)
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
cSTIR_getAcquisitionsStorageScheme()
{
	return charDataHandleFromCharData
		(PETAcquisitionData::storage_scheme().c_str());
}

extern "C"
void* cSTIR_acquisitionsDataFromTemplate(void* ptr_t)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_t, ptr_t);
		shared_ptr<PETAcquisitionData> sptr(sptr_t->new_acquisition_data());
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionsDataFromScannerInfo
(const char* scanner, int span, int max_ring_diff, int view_mash_factor)
{
	std::string storage = PETAcquisitionData::storage_scheme();
	try{
		shared_ptr<ExamInfo> sptr_ei(new ExamInfo());
		if (storage[0] == 'f' || strcmp(storage.c_str(), "default") == 0) {
			shared_ptr<PETAcquisitionDataInFile> 
				sptr(new PETAcquisitionDataInFile
				(sptr_ei, scanner, span, max_ring_diff, view_mash_factor));
			return newObjectHandle(sptr);
		}
		else {
			shared_ptr<PETAcquisitionDataInMemory>
				sptr(new PETAcquisitionDataInMemory
				(sptr_ei, scanner, span, max_ring_diff, view_mash_factor));
			return newObjectHandle(sptr);
		}
	}
	CATCH;
}

extern "C"
void* cSTIR_getAcquisitionsDimensions(const void* ptr_acq, size_t ptr_dim)
{
	try {
		int* dim = (int*)ptr_dim;
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		dim[0] = sptr_ad->get_num_tangential_poss();
		dim[1] = sptr_ad->get_num_views();
		dim[2] = sptr_ad->get_num_sinograms();
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getAcquisitionsData(const void* ptr_acq, size_t ptr_data)
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
void* cSTIR_fillAcquisitionsData(void* ptr_acq, float v)
{
	try {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		sptr_ad->fill(v);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_fillAcquisitionsDataFromAcquisitionsData
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
void* cSTIR_setAcquisitionsData(void* ptr_acq, size_t ptr_data)
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
void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i)
{
	try {
		DataHandle* handle = new DataHandle;
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData* ptr_id = new PETImageData(s);
		shared_ptr<PETImageData> sptr(ptr_id);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		PETImageData* ptr_id = new PETImageData(image);
		shared_ptr<PETImageData> sptr(ptr_id);
		Image3DF& grad = sptr->data();
		if (subset >= 0)
			fun.compute_sub_gradient(grad, image, subset);
		else {
			int nsub = fun.get_num_subsets();
			grad.fill(0.0);
			PETImageData* ptr_id = new PETImageData(image);
			shared_ptr<PETImageData> sptr_sub(ptr_id);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<PETImageData> sptr(new PETImageData(image));
		Image3DF& grad = sptr->data();
		fun.compute_sub_gradient_without_penalty_plus_sensitivity
			(grad, image, subset);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_priorGradient(void* ptr_p, void* ptr_i)
{
	try {
		Prior3DF& prior = objectFromHandle<Prior3DF>(ptr_p);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<PETImageData> sptr(new PETImageData(image));
		Image3DF& grad = sptr->data();
		prior.compute_gradient(grad, image);
		return newObjectHandle(sptr);
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
		shared_ptr<PETImageData> sptr(new PETImageData(voxels));
		return (void*)newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromImage(void* ptr_i)
{
	try {
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<PETImageData> sptr(new PETImageData(image));
		return (void*)newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_writeImage(void* ptr_i, const char* filename)
{
	try {
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		shared_ptr<OutputFileFormat<Image3DF> > format_sptr =
			OutputFileFormat<Image3DF>::default_sptr();
		format_sptr->write_to_file(filename, image);
		return (void*) new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromAcquisitionData(void* ptr_ad)
{
	try {
		shared_ptr<PETAcquisitionData>& sptr_ad =
			objectSptrFromHandle<PETAcquisitionData>(ptr_ad);
		//shared_ptr<ProjDataInfo> sptr_adi =
		//	sptr_ad->get_proj_data_info_sptr();
		shared_ptr<PETImageData> sptr(new PETImageData(*sptr_ad));
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_addShape(void* ptr_i, void* ptr_s, float v)
{
	try {
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
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
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		id.get_voxel_sizes(vs);
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getImageData(const void* ptr_im, size_t ptr_data)
{
	try {
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		float* data = (float*)ptr_data;
		if (id.get_data(data)) {
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
void* cSTIR_setImageData(const void* ptr_im, size_t ptr_data)
{
	try {
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		float* data = (float*)ptr_data;
		if (id.set_data(data)) {
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
void*
cSTIR_norm(const void* ptr_x)
{
	try {
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		return dataHandle(x.norm());
	}
	CATCH;
}

extern "C"
void*
cSTIR_dot(const void* ptr_x, const void* ptr_y)
{
	try {
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		aDataContainer<float>& y =
			objectFromHandle<aDataContainer<float> >(ptr_y);
		return dataHandle(x.dot(y));
	}
	CATCH;
}

extern "C"
void*
cSTIR_mult(float a, const void* ptr_x)
{
	try {
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		shared_ptr<aDataContainer<float> > sptr_z(x.new_data_container());
		//sptr_z->mult(a, x);
		sptr_z->axpby(a, x, 0, x);
		return newObjectHandle<aDataContainer<float> >(sptr_z);
	}
	CATCH;
}

extern "C"
void*
cSTIR_axpby(
	float a, const void* ptr_x,
	float b, const void* ptr_y
) {
	try {
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		aDataContainer<float>& y =
			objectFromHandle<aDataContainer<float> >(ptr_y);
		shared_ptr<aDataContainer<float> > sptr_z(x.new_data_container());
		sptr_z->axpby(a, x, b, y);
		return newObjectHandle<aDataContainer<float> >(sptr_z);
	}
	CATCH;
}

