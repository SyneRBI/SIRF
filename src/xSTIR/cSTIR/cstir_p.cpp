/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

#include <stdio.h>
#include <stdlib.h>

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/STIR/stir_types.h"
#include "sirf/STIR/cstir_p.h"
#include "sirf/STIR/stir_x.h"

using namespace stir;
using namespace sirf;

#define SPTR_FROM_HANDLE(Object, X, H) \
	shared_ptr<Object> X; getObjectSptrFromHandle<Object>(H, X);
#define DYNAMIC_CAST(T, X, Y) T& X = dynamic_cast<T&>(Y)

extern "C"
char* charDataFromHandle(const void* ptr);

static void*
parameterNotFound(const char* name, const char* file, int line) 
{
	DataHandle* handle = new DataHandle;
	std::string error = "parameter ";
	error += name;
	error += " not found";
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

static void*
wrongParameterValue
(const char* name, const char* value, const char* file, int line)
{
	DataHandle* handle = new DataHandle;
	std::string error = "parameter ";
	error += name;
	error += " cannot be ";
	//error += " has wrong value ";
	error += value;
	ExecutionStatus status(error.c_str(), file, line);
	handle->set(0, &status);
	return (void*)handle;
}

static void*
wrongIntParameterValue
(const char* name, int value, const char* file, int line)
{
	char buff[32];
	sprintf(buff, "%d", value);
	return wrongParameterValue(name, buff, file, line);
}

static void*
wrongFloatParameterValue
(const char* name, float value, const char* file, int line)
{
	char buff[32];
	sprintf(buff, "%f", value);
	return wrongParameterValue(name, buff, file, line);
}

void*
sirf::cSTIR_setListmodeToSinogramsParameter(void* hp, const char* name, const void* hv)
{
	ListmodeToSinograms& lm2s = objectFromHandle<ListmodeToSinograms>(hp);
	if (boost::iequals(name, "input"))
		lm2s.set_input(charDataFromHandle(hv));
	else if (boost::iequals(name, "output"))
		lm2s.set_output(charDataFromHandle(hv));
	else if (boost::iequals(name, "template"))
		lm2s.set_template(charDataFromHandle(hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setSeparableGaussianImageFilterParameter
(void* hp, const char* name, const void* hv)
{
	xSTIR_SeparableGaussianImageFilter& f = 
		objectFromHandle<xSTIR_SeparableGaussianImageFilter>(hp);
	if (boost::iequals(name, "fwhms_x")) {
		float v = dataFromHandle<float>(hv);
		f.set_fwhms_xyz('x', v);
	}
	else if (boost::iequals(name, "fwhms_y")) {
		float v = dataFromHandle<float>(hv);
		f.set_fwhms_xyz('y', v);
	}
	else if (boost::iequals(name, "fwhms_z")) {
		float v = dataFromHandle<float>(hv);
		f.set_fwhms_xyz('z', v);
	}
	else if (boost::iequals(name, "max_kernel_size_x")) {
		int mks = dataFromHandle<int>(hv);
		f.set_max_kernel_sizes_xyz('x', mks);
	}
	else if (boost::iequals(name, "max_kernel_size_y")) {
		int mks = dataFromHandle<int>(hv);
		f.set_max_kernel_sizes_xyz('y', mks);
	}
	else if (boost::iequals(name, "max_kernel_size_z")) {
		int mks = dataFromHandle<int>(hv);
		f.set_max_kernel_sizes_xyz('z', mks);
	}
	else if (boost::iequals(name, "normalise")) {
		bool normalise = dataFromHandle<int>(hv);
		f.set_normalise(normalise);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setShapeParameter(void* hp, const char* name, const void* hv)
{
	Shape3D& s = objectFromHandle<Shape3D>(hp);
	Coord3DF origin = s.get_origin();
	float value = dataFromHandle<float>(hv);
	if (boost::iequals(name, "x"))
		origin.x() = value;
	else if (boost::iequals(name, "y"))
		origin.y() = value;
	else if (boost::iequals(name, "z"))
		origin.z() = value;
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	s.set_origin(origin);
	return new DataHandle;
}

void*
sirf::cSTIR_shapeParameter(const DataHandle* handle, const char* name)
{
	Shape3D& s = objectFromHandle<Shape3D>(handle);
	Coord3DF origin = s.get_origin();
	if (boost::iequals(name, "x"))
		return dataHandle<float>(origin.x());
	else if (boost::iequals(name, "y"))
		return dataHandle<float>(origin.y());
	else if (boost::iequals(name, "z"))
		return dataHandle<float>(origin.z());
	else
		return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setEllipsoidalCylinderParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	EllipsoidalCylinder& c =
		objectFromHandle<EllipsoidalCylinder>(hp);
	float value = dataFromHandle<float>(hv);
	if (boost::iequals(name, "length"))
		c.set_length(value);
	else if (boost::iequals(name, "radius_x"))
		c.set_radius_x(value);
	else if (boost::iequals(name, "radius_y"))
		c.set_radius_y(value);
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_ellipsoidalCylinderParameter(const DataHandle* handle, const char* name)
{
	EllipsoidalCylinder& c =
		objectFromHandle<EllipsoidalCylinder>(handle);
	if (boost::iequals(name, "length"))
		return dataHandle<float>(c.get_length());
	if (boost::iequals(name, "radius_x"))
		return dataHandle<float>(c.get_radius_x());
	if (boost::iequals(name, "radius_y"))
		return dataHandle<float>(c.get_radius_y());
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setTruncateToCylindricalFOVImageProcessorParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	CylindricFilter3DF& filter =
		objectFromHandle<CylindricFilter3DF>(hp);
	if (boost::iequals(name, "strictly_less_than_radius"))
		filter.set_strictly_less_than_radius
			(boost::iequals(charDataFromDataHandle(hv), "true"));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_truncateToCylindricalFOVImageProcessorParameter
(const DataHandle* handle, const char* name)
{
	CylindricFilter3DF& filter =
		objectFromHandle<CylindricFilter3DF>(handle);
	if (boost::iequals(name, "strictly_less_than_radius"))
		return dataHandle<int>(filter.get_strictly_less_than_radius());
	else
		return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setRayTracingMatrixParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	RayTracingMatrix& matrix = 
		objectFromHandle<RayTracingMatrix>(hp);
	int value = dataFromHandle<int>(hv);
	if (boost::iequals(name, "num_tangential_LORs"))
		matrix.set_num_tangential_LORs(value);
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_rayTracingMatrixParameter(const DataHandle* handle, const char* name)
{
	RayTracingMatrix& matrix = 
		objectFromHandle<RayTracingMatrix>(handle);
	if (boost::iequals(name, "num_tangential_LORs"))
		return dataHandle<int>(matrix.get_num_tangential_LORs());
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setAcquisitionModelParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	AcqMod3DF& am = objectFromHandle< AcqMod3DF >(hp);
	if (boost::iequals(name, "additive_term")) {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, hv);
		am.set_additive_term(sptr_ad);
	}
	else if (boost::iequals(name, "background_term")) {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, hv);
		am.set_background_term(sptr_ad);
	}
	else if (boost::iequals(name, "asm")) {
		SPTR_FROM_HANDLE(PETAcquisitionSensitivityModel, sptr_asm, hv);
		am.set_asm(sptr_asm);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setAcqModUsingMatrixParameter
(DataHandle* hm, const char* name, const DataHandle* hv)
{
	AcqModUsingMatrix3DF& am = objectFromHandle<AcqModUsingMatrix3DF>(hm);
	if (boost::iequals(name, "matrix")) {
		SPTR_FROM_HANDLE(ProjMatrixByBin, sptr_m, hv);
		am.set_matrix(sptr_m);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_acqModUsingMatrixParameter
(DataHandle* hm, const char* name)
{
	AcqModUsingMatrix3DF& am = objectFromHandle<AcqModUsingMatrix3DF>(hm);
	if (boost::iequals(name, "matrix"))
		return newObjectHandle(am.matrix_sptr());
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setGeneralisedPriorParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	Prior3DF& prior = objectFromHandle< Prior3DF >(hp);
	if (boost::iequals(name, "penalisation_factor"))
		prior.set_penalisation_factor(dataFromHandle<float>((void*)hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_generalisedPriorParameter(const DataHandle* handle, const char* name)
{
	Prior3DF& prior = objectFromHandle< Prior3DF >(handle);
	if (boost::iequals(name, "penalisation_factor"))
		return dataHandle<float>(prior.get_penalisation_factor());
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setQuadraticPriorParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	xSTIR_QuadraticPrior3DF& prior =
		objectFromHandle<xSTIR_QuadraticPrior3DF>(hp);
	if (boost::iequals(name, "only_2D"))
		prior.only2D(dataFromHandle<int>((void*)hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setPLSPriorParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	//xSTIR_PLSPrior3DF& prior =
	//	objectFromHandle<xSTIR_PLSPrior3DF>(hp);
	PLSPrior<float>& prior =
		objectFromHandle<PLSPrior<float> >(hp);
	if (boost::iequals(name, "only_2D"))
		prior.set_only_2D(dataFromHandle<int>((void*)hv));
	else if (boost::iequals(name, "alpha"))
		prior.set_alpha(dataFromHandle<float>((void*)hv));
	else if (boost::iequals(name, "eta"))
		prior.set_eta(dataFromHandle<float>((void*)hv));
	else if (boost::iequals(name, "anatomical_image")) {
		STIRImageData& id = objectFromHandle<STIRImageData>(hv);
		prior.set_anatomical_image_sptr(id.data_sptr());
	}
	else if (boost::iequals(name, "kappa")) {
		STIRImageData& id = objectFromHandle<STIRImageData>(hv);
		prior.set_kappa_sptr(id.data_sptr());
	}
	else if (boost::iequals(name, "kappa_filename"))
		prior.set_kappa_filename(charDataFromDataHandle(hv));
	else if (boost::iequals(name, "anatomical_filename"))
		prior.set_anatomical_filename(charDataFromDataHandle(hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_PLSPriorParameter
(DataHandle* hp, const char* name)
{
	PLSPrior<float>& prior =
		objectFromHandle<PLSPrior<float> >(hp);
	if (boost::iequals(name, "only_2D"))
		return dataHandle<int>(prior.get_only_2D());
	else if (boost::iequals(name, "alpha"))
		return dataHandle<float>(prior.get_alpha());
	else if (boost::iequals(name, "eta"))
		return dataHandle<float>(prior.get_eta());
	else if (boost::iequals(name, "anatomical_image")) {
		sptrImage3DF sptr_im = prior.get_anatomical_image_sptr();
		shared_ptr<STIRImageData> sptr_id(new STIRImageData(sptr_im));
		return newObjectHandle(sptr_id);
	}
	else if (boost::iequals(name, "kappa")) {
		sptrImage3DF sptr_im = prior.get_kappa_sptr();
		shared_ptr<STIRImageData> sptr_id(new STIRImageData(sptr_im));
		return newObjectHandle(sptr_id);
	}
	else if (boost::iequals(name, "norm")) {
		sptrImage3DF sptr_im = prior.get_norm_sptr();
		shared_ptr<STIRImageData> sptr_id(new STIRImageData(sptr_im));
		return newObjectHandle(sptr_id);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setGeneralisedObjectiveFunctionParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	ObjectiveFunction3DF& obj_fun =
		objectFromHandle< ObjectiveFunction3DF >(hp);
	if (boost::iequals(name, "prior")) {
		SPTR_FROM_HANDLE(GeneralisedPrior<Image3DF>, sptr_p, hv);
		obj_fun.set_prior_sptr(sptr_p);
	}
	else if (boost::iequals(name, "num_subsets"))
		obj_fun.set_num_subsets(dataFromHandle<int>((void*)hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_generalisedObjectiveFunctionParameter
(const DataHandle* handle, const char* name)
{
	ObjectiveFunction3DF& obj_fun =
		objectFromHandle< ObjectiveFunction3DF >(handle);
	if (boost::iequals(name, "prior"))
		return newObjectHandle(obj_fun.get_prior_sptr());
	else if (boost::iequals(name, "num_subsets"))
		return dataHandle<int>(obj_fun.get_num_subsets());
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	PoissonLogLhLinModMean3DF& obj_fun = 
		objectFromHandle<PoissonLogLhLinModMean3DF>(hp);
	if (boost::iequals(name, "sensitivity_filename"))
		obj_fun.set_sensitivity_filename(charDataFromDataHandle(hv));
	else if (boost::iequals(name, "use_subset_sensitivities"))
		obj_fun.set_use_subset_sensitivities
			(boost::iequals(charDataFromDataHandle(hv), "true"));
	else if (boost::iequals(name, "recompute_sensitivity"))
		obj_fun.set_recompute_sensitivity
			(boost::iequals(charDataFromDataHandle(hv), "true"));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setPoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF&
		obj_fun = objectFromHandle
		<xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF>
		(hp);
	if (boost::iequals(name, "input_filename"))
		obj_fun.set_input_file(charDataFromDataHandle(hv));
	else if (boost::iequals(name, "zero_seg0_end_planes"))
		obj_fun.set_zero_seg0_end_planes
			(boost::iequals(charDataFromDataHandle(hv), "true"));
	//else if (boost::iequals(name, "max_segment_num_to_process"))
	//	obj_fun.set_max_segment_num_to_process(dataFromHandle<int>((void*)hv));
	else if (boost::iequals(name, "acquisition_data")) {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, hv);
		obj_fun.set_acquisition_data(sptr_ad);
	}
	else if (boost::iequals(name, "acquisition_model")) {
		SPTR_FROM_HANDLE(AcqMod3DF, sptr_am, hv);
		obj_fun.set_acquisition_model(sptr_am);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjDataParameter
(const DataHandle* handle, const char* name)
{
	xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF&
		obj_fun = objectFromHandle
		<xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF>
		(handle);
	//if (boost::iequals(name, "projector_pair_type"))
	//	return newObjectHandle(obj_fun.get_projector_pair_sptr());
	if (boost::iequals(name, "acquisition_model"))
		return newObjectHandle(obj_fun.acquisition_model_sptr());
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setReconstructionParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	Reconstruction<Image3DF>& recon =
		objectFromHandle< Reconstruction<Image3DF> >(hp);
	if (boost::iequals(name, "output_filename_prefix"))
		recon.set_output_filename_prefix(charDataFromDataHandle(hv));
	else if (boost::iequals(name, "input_data")) {
		SPTR_FROM_HANDLE(PETAcquisitionData, sptr_ad, hv);
		recon.set_input_data(sptr_ad->data()); // objectSptrFromHandle<PETAcquisitionData>(hv)->data());
	}
	else if (boost::iequals(name, "disable_output")) {
		recon.set_disable_output(true);
	}
	else if (boost::iequals(name, "enable_output")) {
		recon.set_disable_output(false);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_setIterativeReconstructionParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	IterativeReconstruction3DF& recon =
		objectFromHandle<IterativeReconstruction3DF>(hp);
	if (boost::iequals(name, "inter_iteration_filter_type")) {
		SPTR_FROM_HANDLE(DataProcessor3DF, sptr_f, hv);
		recon.set_inter_iteration_filter_ptr(sptr_f);
	}
	else if (boost::iequals(name, "objective_function")) {
		SPTR_FROM_HANDLE(ObjectiveFunction3DF, sptr_obf, hv);
		recon.set_objective_function_sptr(sptr_obf);
	}
	else if (boost::iequals(name, "initial_estimate")) {
		//xSTIR_IterativeReconstruction3DF& xrecon =
		//	(xSTIR_IterativeReconstruction3DF&)(recon);
		DYNAMIC_CAST(xSTIR_IterativeReconstruction3DF, xrecon, recon);
		xrecon.set_initial_estimate_file(charDataFromDataHandle(hv));
	}
	else {
		int value = dataFromHandle<int>((void*)hv);
		if (boost::iequals(name, "num_subsets"))
			recon.set_num_subsets(value);
		else if (boost::iequals(name, "start_subset_num"))
			recon.set_start_subset_num(value);
		else if (boost::iequals(name, "num_subiterations"))
			recon.set_num_subiterations(value);
		else if (boost::iequals(name, "start_subiteration_num"))
			recon.set_start_subiteration_num(value);
		else if (boost::iequals(name, "subiteration_num")) {
			//xSTIR_IterativeReconstruction3DF& xrecon =
			//	(xSTIR_IterativeReconstruction3DF&)(recon);
			DYNAMIC_CAST(xSTIR_IterativeReconstruction3DF, xrecon, recon);
			xrecon.subiteration() = value;
		}
		else if (boost::iequals(name, "save_interval"))
			recon.set_save_interval(value);
		else if (boost::iequals(name, "inter_iteration_filter_interval"))
			recon.set_inter_iteration_filter_interval(value);
		else
			return parameterNotFound(name, __FILE__, __LINE__);
	}
	return new DataHandle;
}

void*
sirf::cSTIR_iterativeReconstructionParameter
(const DataHandle* handle, const char* name)
{
	IterativeReconstruction3DF& recon = 
		objectFromHandle<IterativeReconstruction3DF>(handle);
	if (boost::iequals(name, "num_subsets"))
		return dataHandle<int>(recon.get_num_subsets());
	if (boost::iequals(name, "start_subset_num"))
		return dataHandle<int>(recon.get_start_subset_num());
	if (boost::iequals(name, "num_subiterations"))
		return dataHandle<int>(recon.get_num_subiterations());
	if (boost::iequals(name, "start_subiteration_num"))
		return dataHandle<int>(recon.get_start_subiteration_num());
	if (boost::iequals(name, "subiteration_num")) {
		//xSTIR_IterativeReconstruction3DF& xrecon =
		//	(xSTIR_IterativeReconstruction3DF&)(recon);
		//DYNAMIC_CAST(xSTIR_IterativeReconstruction3DF, xrecon, recon);
		int iter = recon.get_subiteration_num();
		return dataHandle<int>(iter);
		//return dataHandle<int>(xrecon.subiteration());
	}
	if (boost::iequals(name, "objective_function"))
		return newObjectHandle(recon.get_objective_function_sptr());
	if (boost::iequals(name, "inter_iteration_filter_type"))
		return newObjectHandle(recon.get_inter_iteration_filter_sptr());
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setOSMAPOSLParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	OSMAPOSLReconstruction<Image3DF>& recon =
		objectFromHandle<OSMAPOSLReconstruction<Image3DF> >(hp);
    if (boost::iequals(name, "set_maximum_relative_change"))
            recon.set_maximum_relative_change(dataFromHandle<double>((void*)hv));
    if (boost::iequals(name, "set_minimum_relative_change"))
            recon.set_minimum_relative_change(dataFromHandle<double>((void*)hv));
	if (boost::iequals(name, "MAP_model"))
		recon.set_MAP_model(charDataFromDataHandle(hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void* 
sirf::cSTIR_OSMAPOSLParameter(const DataHandle* handle, const char* name)
{
	OSMAPOSLReconstruction<Image3DF>& recon =
		objectFromHandle<OSMAPOSLReconstruction<Image3DF> >(handle);
	if (boost::iequals(name, "objective_function"))
		return newObjectHandle(recon.get_objective_function_sptr());
	return parameterNotFound(name, __FILE__, __LINE__);
}
#ifdef USE_HKEM
void*
sirf::cSTIR_setKOSMAPOSLParameter
(DataHandle* hp, const char* name, const DataHandle* hv)
{
	KOSMAPOSLReconstruction<Image3DF>& recon =
		objectFromHandle<KOSMAPOSLReconstruction<Image3DF> >(hp);
	if (boost::iequals(name, "anatomical_prior")) {
		SPTR_FROM_HANDLE(STIRImageData, sptr_ap, hv);
		recon.set_anatomical_prior_sptr(sptr_ap->data_sptr());
	}
	else if (boost::iequals(name, "num_neighbours")) {
		int value = dataFromHandle<int>((void*)hv);
		recon.set_num_neighbours(value);
	}
	else if (boost::iequals(name, "num_non_zero_features"))
		recon.set_num_non_zero_feat(dataFromHandle<int>((void*)hv));
	else if (boost::iequals(name, "sigma_m"))
		recon.set_sigma_m(dataFromHandle<float>((void*)hv));
	else if (boost::iequals(name, "sigma_p"))
		recon.set_sigma_p(dataFromHandle<float>((void*)hv));
	else if (boost::iequals(name, "sigma_dp"))
		recon.set_sigma_dp(dataFromHandle<float>((void*)hv));
	else if (boost::iequals(name, "sigma_dm"))
		recon.set_sigma_dm(dataFromHandle<float>((void*)hv));
	else if (boost::iequals(name, "only_2D"))
		recon.set_only_2D(dataFromHandle<int>((void*)hv));
	else if (boost::iequals(name, "hybrid"))
		recon.set_hybrid(dataFromHandle<int>((void*)hv));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}
#endif
void*
sirf::cSTIR_setOSSPSParameter(DataHandle* hp, const char* name, const DataHandle* hv)
{
	xSTIR_OSSPSReconstruction3DF& recon =
		objectFromHandle<xSTIR_OSSPSReconstruction3DF >(hp);
	if (boost::iequals(name, "relaxation_parameter"))
		recon.relaxation_parameter_value() = dataFromHandle<float>(hv);
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void* 
sirf::cSTIR_OSSPSParameter(const DataHandle* handle, const char* name)
{
	xSTIR_OSSPSReconstruction3DF& recon =
		objectFromHandle<xSTIR_OSSPSReconstruction3DF>(handle);
	return parameterNotFound(name, __FILE__, __LINE__);
}

void*
sirf::cSTIR_setFBP2DParameter(DataHandle* hp, const char* name, const DataHandle* hv)
{
	xSTIR_FBP2DReconstruction& recon =
		objectFromHandle<xSTIR_FBP2DReconstruction >(hp);
	if (boost::iequals(name, "input")) {
		PETAcquisitionData& acq_data = objectFromHandle<PETAcquisitionData>(hv);
		recon.set_input(acq_data);
	}
	else if (boost::iequals(name, "zoom")) {
		double zoom = dataFromHandle<float>(hv);
		recon.set_zoom(zoom);
	}
	else if (boost::iequals(name, "xy")) {
		int xy = dataFromHandle<int>(hv);
		recon.set_output_image_size_xy(xy);
	}
	else if (boost::iequals(name, "alpha")) {
		double alpha = dataFromHandle<float>(hv);
		recon.set_alpha_ramp(alpha);
	}
	else if (boost::iequals(name, "fc")) {
		double fc = dataFromHandle<float>(hv);
		recon.set_frequency_cut_off(fc);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}

void*
sirf::cSTIR_FBP2DParameter(DataHandle* hp, const char* name)
{
	xSTIR_FBP2DReconstruction& recon =
		objectFromHandle<xSTIR_FBP2DReconstruction >(hp);
	if (boost::iequals(name, "output"))
		return newObjectHandle(recon.get_output());
	return parameterNotFound(name, __FILE__, __LINE__);
}

