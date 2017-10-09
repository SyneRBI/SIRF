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
cSTIR_newReconstructionMethod(const char* parFile)
{
	try {
		if (strlen(parFile) > 0) {
			NEW_SPTR(Reconstruction<Image3DF>, ptr_sptr, Method(parFile));
			return newObjectHandle(ptr_sptr);
		}
		else {
			NEW_SPTR(Reconstruction<Image3DF>, ptr_sptr, Method);
			return newObjectHandle(ptr_sptr);
		}
	}
	CATCH;
}

extern "C"
void* cSTIR_newObject(const char* name)
{
	try {
		if (boost::iequals(name,
			"PoissonLogLikelihoodWithLinearModelForMeanAndProjData"))
			return newObjectHandle
			//< ObjectiveFunction3DF, 
			<xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF>();
		if (boost::iequals(name, "AcqModUsingMatrix"))
			return newObjectHandle
			//< AcqMod3DF, 
			<AcqModUsingMatrix3DF>();
		if (boost::iequals(name, "RayTracingMatrix"))
			return newObjectHandle
			//< ProjMatrixByBin, 
			<RayTracingMatrix>();
		if (boost::iequals(name, "QuadraticPrior"))
			return newObjectHandle
			//< Prior3DF, 
			<QuadPrior3DF>();
		if (boost::iequals(name, "TruncateToCylindricalFOVImageProcessor"))
			return newObjectHandle
			//< DataProcessor3DF, 
			<CylindricFilter3DF>();
		if (boost::iequals(name, "EllipsoidalCylinder"))
			return newObjectHandle
			//< Shape3D, 
			<EllipsoidalCylinder>();
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
		if (boost::iequals(obj, "Shape"))
			return cSTIR_setShapeParameter(hs, name, hv);
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
			//sptrImage3DF* ptr_sptr = new sptrImage3DF
			//	(read_from_file<Image3DF>(filename));
			PETImageData* ptr_id = 
				new PETImageData(read_from_file<Image3DF>(filename));
			shared_ptr<PETImageData>* ptr_sptr =
				new shared_ptr<PETImageData>(ptr_id);
			return newObjectHandle(ptr_sptr);
		}
		if (boost::iequals(name, "AcquisitionData")) {
			//writeText("\nreading ");
			//writeText(filename);
			//NEW(shared_ptr<ProjData>, ptr_sptr);
			//*ptr_sptr = boost::static_pointer_cast<PETAcquisitionData>
			//	(ProjData::read_from_file(filename));
			//writeText("ok\n");
			NEW_SPTR(PETAcquisitionData, ptr_sptr, 
				PETAcquisitionDataInFile(filename));
			//if (ptr_sptr->get())
			//	std::cout << "ok\n";
			//else
			//	std::cout << "read_from_file failed\n";
			return newObjectHandle(ptr_sptr);
		}
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

//extern "C"
//void* cSTIR_setupObject(const char* obj, void* ptr_obj)
//{
//	try {
//		bool status = 1;
//		if (boost::iequals(obj, "GeneralisedObjectiveFunction")) {
//			xSTIR_GeneralisedObjectiveFunction3DF& obj_fun =
//				objectFromHandle<xSTIR_GeneralisedObjectiveFunction3DF>(ptr_obj);
//			status = obj_fun.post_process();
//		}
//		//else if (boost::iequals(obj, "GeneralisedPrior"))
//		//		status = xSTIR_setupPrior(ho->data());
//		DataHandle* handle = new DataHandle;
//		if (status) {
//			ExecutionStatus status
//				("cSTIR_setupObject failed", __FILE__, __LINE__);
//			handle->set(0, &status);
//		}
//		return (void*)handle;
//	}
//	CATCH;
//}

extern "C"
void* cSTIR_applyDataProcessor(const void* ptr_p, void* ptr_i)
{
	try {
		DataProcessor<Image3DF>& processor =
			objectFromHandle<DataProcessor<Image3DF> >(ptr_p);
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		processor.apply(image);
		return (void*) new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupAcquisitionModel(void* ptr_am, void* ptr_dt, void* ptr_im)
{
	try {
		//writeText("setting up acquisition model\n");
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		//sptrProjData sptr_dt = objectSptrFromHandle<ProjData>(ptr_dt);
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_dt, ptr_dt);
		//sptrImage3DF sptr_im = objectSptrFromHandle<Image3DF>(ptr_im);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		sptrImage3DF sptr_im = id.data_sptr();
		//std::cout << "setting up acquisition model...\n";
		Succeeded s = am.set_up(sptr_dt, sptr_im);
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
		//Image3DF& im = objectFromHandle<Image3DF>(ptr_im);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		Image3DF& im = id.data();
		//sptrProjData* ptr_sptr = new sptrProjData;
		//*ptr_sptr = am.forward(im, datafile);
		//return newObjectHandle(ptr_sptr);
		NEW(shared_ptr<PETAcquisitionData>, sptr_ad);
		*sptr_ad = am.forward(im);
		return newObjectHandle(sptr_ad);
	}
	CATCH;
}

extern "C"
void* cSTIR_acquisitionModelBwd(void* ptr_am, void* ptr_ad)
{
	try {
		AcqMod3DF& am = objectFromHandle<AcqMod3DF>(ptr_am);
		//ProjData& ad = objectFromHandle<ProjData>(ptr_ad);
		PETAcquisitionData& ad = objectFromHandle<PETAcquisitionData>(ptr_ad);
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(am.backward(*ad.data()));
		PETImageData* ptr_id = new PETImageData(am.backward(*ad.data()));
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		return newObjectHandle(ptr_sptr);
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
	//PETAcquisitionData::set_storage_scheme(scheme);
	//return (void*)new DataHandle;
}

extern "C"
void* cSTIR_acquisitionsDataFromTemplate(void* ptr_t)
{
	try {
		//sptrProjData& sptr_t = objectSptrFromHandle<ProjData>(ptr_t);
		//NEW_SPTR(ProjData, ptr_sptr,
		//	ProjDataInMemory(sptr_t->get_exam_info_sptr(),
		//					 sptr_t->get_proj_data_info_sptr()));
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_t, ptr_t);
		NEW_SPTR_FROM_PTR
			(PETAcquisitionData, ptr_sptr, sptr_t->new_acquisition_data());
		return newObjectHandle(ptr_sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_getAcquisitionsDimensions(const void* ptr_acq, size_t ptr_dim)
{
	try {
		int* dim = (int*)ptr_dim;
		//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(ptr_acq);
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
		//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(ptr_acq);
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
		//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(ptr_acq);
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
		//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(ptr_acq);
		//sptrProjData& sptr_from = objectSptrFromHandle<ProjData>(ptr_from);
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
		//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(ptr_acq);
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, ptr_acq);
		float *data = (float *)ptr_data;
		sptr_ad->fill_from(data);
		return (void*)new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setupReconstruction(void* ptr_r, void* ptr_i)
{
	try {
		DataHandle* handle = new DataHandle;
		//sptrImage3DF& sptr_image = objectSptrFromHandle<Image3DF>(ptr_i);
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
		//sptrImage3DF& sptr_image = objectSptrFromHandle<Image3DF>(ptr_i);
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
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
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
		//sptrImage3DF& sptr_image = objectSptrFromHandle<Image3DF>(ptr_i);
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
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
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
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(s.clone());
		PETImageData* ptr_id = new PETImageData(s);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		return newObjectHandle(ptr_sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_objectiveFunctionGradient(void* ptr_f, void* ptr_i, int subset)
{
	try {
		ObjectiveFunction3DF& fun = objectFromHandle< ObjectiveFunction3DF>(ptr_f);
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(image.get_empty_copy());
		PETImageData* ptr_id = new PETImageData(image);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		//Image3DF& grad = **ptr_sptr;
		Image3DF& grad = (**ptr_sptr).data();
		if (subset >= 0)
			fun.compute_sub_gradient(grad, image, subset);
		else {
			int nsub = fun.get_num_subsets();
			grad.fill(0.0);
			//sptrImage3DF* sptr_subgrad = new sptrImage3DF(image.get_empty_copy());
			//Image3DF& subgrad = **sptr_subgrad;
			PETImageData* ptr_id = new PETImageData(image);
			NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
			Image3DF& subgrad = (**ptr_sptr).data();
			for (int sub = 0; sub < nsub; sub++) {
				fun.compute_sub_gradient(subgrad, image, sub);
				grad += subgrad;
			}
		}
		return newObjectHandle(ptr_sptr);
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
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(image.get_empty_copy());
		//Image3DF& grad = **ptr_sptr; // ->get();
		PETImageData* ptr_id = new PETImageData(image);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		Image3DF& grad = (**ptr_sptr).data();
		fun.compute_sub_gradient_without_penalty_plus_sensitivity
			(grad, image, subset);
		return newObjectHandle(ptr_sptr);
	}
	CATCH;
}

extern "C"
void*
cSTIR_priorGradient(void* ptr_p, void* ptr_i)
{
	try {
		Prior3DF& prior = objectFromHandle<Prior3DF>(ptr_p);
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(image.get_empty_copy());
		//Image3DF& grad = **ptr_sptr; // ->get();
		PETImageData* ptr_id = new PETImageData(image);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		Image3DF& grad = (**ptr_sptr).data();
		prior.compute_gradient(grad, image);
		return newObjectHandle(ptr_sptr);
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
		sptrVoxels3DF* sptr = new sptrVoxels3DF(
			new Voxels3DF(IndexRange3D(0, nz - 1,
			-(ny / 2), -(ny / 2) + ny - 1, -(nx / 2), -(nx / 2) + nx - 1),
			Coord3DF(z, y, x),
			Coord3DF(sz, sy, sx)));
		(*sptr)->fill(0);
		return newObjectHandle(sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromVoxels(void* ptr_v)
{
	try {
		Voxels3DF& voxels = objectFromHandle<Voxels3DF>(ptr_v);
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(voxels.clone());
		PETImageData* ptr_id = new PETImageData(voxels);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		return newObjectHandle(ptr_sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_imageFromImage(void* ptr_i)
{
	try {
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
		PETImageData& id = objectFromHandle<PETImageData>(ptr_i);
		Image3DF& image = id.data();
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(image.clone());
		PETImageData* ptr_id = new PETImageData(image);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		return newObjectHandle(ptr_sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_writeImage(void* ptr_i, const char* filename)
{
	try {
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
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
		//sptrProjData& sptr_ad = objectSptrFromHandle<ProjData>(ptr_ad);
		shared_ptr<PETAcquisitionData>& sptr_ad =
			objectSptrFromHandle<PETAcquisitionData>(ptr_ad);
		shared_ptr<ProjDataInfo> sptr_adi =
			sptr_ad->get_proj_data_info_sptr();
		//Voxels3DF* ptr_voxels = new Voxels3DF(*sptr_adi);
		//sptrImage3DF* ptr_sptr = new sptrImage3DF(ptr_voxels);
		PETImageData* ptr_id = new PETImageData(*sptr_adi);
		NEW_SPTR_FROM_PTR(PETImageData, ptr_sptr, ptr_id);
		return newObjectHandle(ptr_sptr);
	}
	CATCH;
}

extern "C"
void* cSTIR_addShape(void* ptr_i, void* ptr_s, float v)
{
	try {
		//Image3DF& image = objectFromHandle<Image3DF>(ptr_i);
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
		//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>(ptr_i);
		//if (ptr_image == 0)
		//	return;
		//Image3DF& image = *ptr_image;
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
		dim[0] = 0;
		dim[1] = 0;
		dim[2] = 0;
		//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>(ptr_im);
		//if (ptr_image == 0)
		//	return;
		//Image3DF& image = *ptr_image;
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		Image3DF& image = id.data();
		Coordinate3D<int> min_indices;
		Coordinate3D<int> max_indices;
		if (!image.get_regular_range(min_indices, max_indices)) {
			ExecutionStatus status("not a regular image", __FILE__, __LINE__);
			DataHandle* handle = new DataHandle;
			handle->set(0, &status);
			return (void*)handle;
		}
		image.get_regular_range(min_indices, max_indices);
		for (int i = 0; i < 3; i++)
			dim[i] = max_indices[i + 1] - min_indices[i + 1] + 1;
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_getImageData(const void* ptr_im, size_t ptr_data) 
{
	try {
		//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>(ptr_im);
		//if (ptr_image == 0)
		//	return;
		//Image3DF& image = *ptr_image;
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		Image3DF& image = id.data();
		Coordinate3D<int> min_indices;
		Coordinate3D<int> max_indices;
		float* data = (float*)ptr_data;
		if (!image.get_regular_range(min_indices, max_indices)) {
			ExecutionStatus status("not a regular image", __FILE__, __LINE__);
			DataHandle* handle = new DataHandle;
			handle->set(0, &status);
			return (void*)handle;
		}
		for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
			for (int y = min_indices[2]; y <= max_indices[2]; y++) {
				for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
					data[i] = image[z][y][x];
				}
			}
		}
		return new DataHandle;
	}
	CATCH;
}

extern "C"
void* cSTIR_setImageData(const void* ptr_im, size_t ptr_data) 
{
	try {
		//Image3DF* ptr_image = objectPtrFromHandle<Image3DF>(ptr_im);
		//if (ptr_image == 0)
		//	return;
		//Image3DF& image = *ptr_image;
		PETImageData& id = objectFromHandle<PETImageData>(ptr_im);
		Image3DF& image = id.data();
		Coordinate3D<int> min_indices;
		Coordinate3D<int> max_indices;
		float* data = (float*)ptr_data;
		if (!image.get_regular_range(min_indices, max_indices)) {
			ExecutionStatus status("not a regular image", __FILE__, __LINE__);
			DataHandle* handle = new DataHandle;
			handle->set(0, &status);
			return (void*)handle;
		}
		for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
			for (int y = min_indices[2]; y <= max_indices[2]; y++) {
				for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
					image[z][y][x] = data[i];
				}
			}
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
		//CAST_PTR(DataHandle, h_x, ptr_x);
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		//objectFromHandle<aDataContainer<float> >(h_x);
		float* result = (float*)malloc(sizeof(float));
		*result = x.norm();
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cSTIR_dot(const void* ptr_x, const void* ptr_y)
{
	try {
		//CAST_PTR(DataHandle, h_x, ptr_x);
		//CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		aDataContainer<float>& y =
			objectFromHandle<aDataContainer<float> >(ptr_y);
		float* result = (float*)malloc(sizeof(float));
		*result = x.dot(y);
		DataHandle* handle = new DataHandle;
		handle->set(result, 0, GRAB);
		return (void*)handle;
	}
	CATCH;
}

extern "C"
void*
cSTIR_mult(float a, const void* ptr_x)
{
	try {
		//CAST_PTR(DataHandle, h_x, ptr_x);
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		//shared_ptr<aDataContainer<float> > sptr_z =
		//	x.new_data_container();
		shared_ptr<aDataContainer<float> > sptr_z(x.new_data_container());
		//sptr_z->mult(a, x);
		sptr_z->axpby(a, x, 0, x);
		return sptrObjectHandle<aDataContainer<float> >(sptr_z);
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
		//CAST_PTR(DataHandle, h_x, ptr_x);
		//CAST_PTR(DataHandle, h_y, ptr_y);
		aDataContainer<float>& x =
			objectFromHandle<aDataContainer<float> >(ptr_x);
		aDataContainer<float>& y =
			objectFromHandle<aDataContainer<float> >(ptr_y);
		shared_ptr<aDataContainer<float> > sptr_z(x.new_data_container());
		//shared_ptr<aDataContainer<float> > sptr_z =
		//	x.new_data_container();
		sptr_z->axpby(a, x, b, y);
		return sptrObjectHandle<aDataContainer<float> >(sptr_z);
	}
	CATCH;
}

