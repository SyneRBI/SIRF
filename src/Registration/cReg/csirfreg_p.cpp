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

#include "stir_x.h"
#include "data_handle.h"
#include "csirfreg_p.h"
#include "SIRFImageData.h"
#include "SIRFReg.h"
#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFRegImageWeightedMean.h"

using namespace stir;
using namespace sirf;

extern "C"
char* charDataFromHandle(const void* ptr);
extern "C"
int intDataFromHandle(const void* ptr);

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

// ------------------------------------------------------------------------------------ //
//   SIRFImageData
// ------------------------------------------------------------------------------------ //
void*
sirf::cSIRFReg_SIRFImageDataParameter(const DataHandle* handle, const char* name)
{
	SIRFImageData& s = objectFromHandle<SIRFImageData>(handle);
	if (boost::iequals(name, "max"))
		return dataHandle<float>(s.get_max());
	if (boost::iequals(name, "min"))
		return dataHandle<float>(s.get_min());
	else
		return parameterNotFound(name, __FILE__, __LINE__);
}
// ------------------------------------------------------------------------------------ //
//   SIRFReg
// ------------------------------------------------------------------------------------ //
// set
void*
sirf::cSIRFReg_setSIRFRegParameter(void* hp, const char* name, const void* hv)
{
	SIRFReg& s = objectFromHandle<SIRFReg>(hp);
	if (boost::iequals(name, "parameter_file"))
		s.set_parameter_file(charDataFromHandle(hv));
	else if (boost::iequals(name, "reference_image")) {
		const SIRFImageData& im = objectFromHandle<const SIRFImageData>(hv);
		s.set_reference_image(im);
	}
	else if (boost::iequals(name, "floating_image")) {
		const SIRFImageData& im = objectFromHandle<const SIRFImageData>(hv);
		s.set_floating_image(im);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}
// get
void*
sirf::cSIRFReg_SIRFRegParameter(const DataHandle* handle, const char* name)
{
	SIRFReg& s = objectFromHandle<SIRFReg>(handle);
	if (boost::iequals(name, "output")) {
        shared_ptr<SIRFImageData> sptr_id(new SIRFImageData(s.get_output()));
        return newObjectHandle(sptr_id);
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
}
// ------------------------------------------------------------------------------------ //
//   SIRFRegNiftyF3dSym
// ------------------------------------------------------------------------------------ //
// set
void*
sirf::cSIRFReg_setSIRFRegNiftyF3dSymParameter(void* hp, const char* name, const void* hv)
{
    SIRFRegNiftyF3dSym<float>& s = objectFromHandle<SIRFRegNiftyF3dSym<float> >(hp);
    if (boost::iequals(name, "floating_time_point"))
        s.set_floating_time_point(intDataFromHandle(hv));
    else if (boost::iequals(name, "reference_time_point"))
        s.set_reference_time_point(intDataFromHandle(hv));
    else if (boost::iequals(name, "initial_affine_transformation"))
        s.set_initial_affine_transformation(charDataFromHandle(hv));
    else
        return parameterNotFound(name, __FILE__, __LINE__);
    return new DataHandle;
}
// ------------------------------------------------------------------------------------ //
//   SIRFRegNiftyResample
// ------------------------------------------------------------------------------------ //
// set
void*
sirf::cSIRFReg_setSIRFRegNiftyResampleParameter(void* hp, const char* name, const void* hv)
{
    SIRFRegNiftyResample& s = objectFromHandle<SIRFRegNiftyResample>(hp);
    if (boost::iequals(name, "reference_image"))
        s.set_reference_image(objectFromHandle<SIRFImageData>(hv));
    else if (boost::iequals(name, "floating_image"))
        s.set_floating_image(objectFromHandle<SIRFImageData>(hv));
    else if (boost::iequals(name, "transformation_matrix"))
        s.set_transformation_matrix(charDataFromHandle(hv));
    else if (boost::iequals(name, "displacement_field"))
        s.set_displacement_field(objectFromHandle<SIRFImageDataDeformation>(hv));
    else if (boost::iequals(name, "deformation_field"))
        s.set_deformation_field(objectFromHandle<SIRFImageDataDeformation>(hv));
    else if (boost::iequals(name, "interpolation_type"))
        s.set_interpolation_type(intDataFromHandle(hv));
    else
        return parameterNotFound(name, __FILE__, __LINE__);
    return new DataHandle;
}
// get
void*
sirf::cSIRFReg_SIRFRegNiftyResampleParameter(const DataHandle* handle, const char* name)
{
    SIRFRegNiftyResample& s = objectFromHandle<SIRFRegNiftyResample>(handle);
    if (boost::iequals(name, "output")) {
        shared_ptr<SIRFImageData> sptr_id(new SIRFImageData(s.get_output()));
        return newObjectHandle(sptr_id);
    }
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}

// ------------------------------------------------------------------------------------ //
//   SIRFRegImageWeightedMean3D
// ------------------------------------------------------------------------------------ //
// get
void*
sirf::cSIRFReg_SIRFRegImageWeightedMean3DParameter(const DataHandle* handle, const char* name)
{
    SIRFRegImageWeightedMean3D& s = objectFromHandle<SIRFRegImageWeightedMean3D>(handle);
    if (boost::iequals(name, "output")) {
        shared_ptr<SIRFImageData> sptr_id(new SIRFImageData(s.get_output()));
        return newObjectHandle(sptr_id);
    }
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}


// ------------------------------------------------------------------------------------ //
//   SIRFRegImageWeightedMean3D
// ------------------------------------------------------------------------------------ //
// get
void*
sirf::cSIRFReg_SIRFRegImageWeightedMean4DParameter(const DataHandle* handle, const char* name)
{
    SIRFRegImageWeightedMean4D& s = objectFromHandle<SIRFRegImageWeightedMean4D>(handle);
    if (boost::iequals(name, "output")) {
        shared_ptr<SIRFImageDataDeformation> sptr_id(new SIRFImageDataDeformation(s.get_output()));
        return newObjectHandle(sptr_id);
    }
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}

