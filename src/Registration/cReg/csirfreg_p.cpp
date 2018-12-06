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

#include "data_handle.h"
#include "csirfreg_p.h"
#include "NiftiImageData3D.h"
#include "SIRFReg.h"
#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFRegImageWeightedMean.h"

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
//   NiftiImageData
// ------------------------------------------------------------------------------------ //
void*
sirf::cSIRFReg_NiftiImageDataParameter(const DataHandle* handle, const char* name)
{
    NiftiImageData<float>& s = objectFromHandle<NiftiImageData<float> >(handle);
    if (strcmp(name, "max") == 0)
        return dataHandle<float>(s.get_max());
    if (strcmp(name, "min") == 0)
        return dataHandle<float>(s.get_min());
    if (strcmp(name, "sum") == 0)
        return dataHandle<float>(s.get_sum());
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
    SIRFReg<float>& s = objectFromHandle<SIRFReg<float> >(hp);
    if (strcmp(name, "parameter_file") == 0)
        s.set_parameter_file(charDataFromHandle(hv));
    else if (strcmp(name, "reference_image") == 0)
        s.set_reference_image(std::make_shared<const NiftiImageData3D<float> >(objectFromHandle<const NiftiImageData3D<float> >(hv)));
    else if (strcmp(name, "floating_image") == 0)
        s.set_floating_image(std::make_shared<const NiftiImageData3D<float> >(objectFromHandle<const NiftiImageData3D<float> >(hv)));
    else if (strcmp(name, "reference_mask") == 0)
        s.set_reference_mask(std::make_shared<const NiftiImageData3D<float> >(objectFromHandle<const NiftiImageData3D<float> >(hv)));
    else if (strcmp(name, "floating_mask") == 0)
        s.set_floating_mask(std::make_shared<const NiftiImageData3D<float> >(objectFromHandle<const NiftiImageData3D<float> >(hv)));
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}
// get
void*
sirf::cSIRFReg_SIRFRegParameter(const DataHandle* handle, const char* name)
{
    SIRFReg<float>& s = objectFromHandle<SIRFReg<float> >(handle);
    if (strcmp(name, "output") == 0) {
        return newObjectHandle(s.get_output());
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
    if (strcmp(name, "floating_time_point") == 0)
        s.set_floating_time_point(intDataFromHandle(hv));
    else if (strcmp(name, "reference_time_point") == 0)
        s.set_reference_time_point(intDataFromHandle(hv));
    else if (strcmp(name, "initial_affine_transformation") == 0) {
        const std::shared_ptr<const SIRFRegAffineTransformation<float> > mat_sptr =
                std::make_shared<const SIRFRegAffineTransformation<float> >(objectFromHandle<const SIRFRegAffineTransformation<float> >(hv));
        s.set_initial_affine_transformation(mat_sptr);
    }
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
    SIRFRegNiftyResample<float>& s = objectFromHandle<SIRFRegNiftyResample<float> >(hp);
    if (strcmp(name, "reference_image") == 0)
        s.set_reference_image(std::make_shared<const NiftiImageData3D<float> >(objectFromHandle<const NiftiImageData3D<float> >(hv)));
    else if (strcmp(name, "floating_image") == 0)
        s.set_floating_image(std::make_shared<const NiftiImageData3D<float> >(objectFromHandle<const NiftiImageData3D<float> >(hv)));
    else if (strcmp(name, "interpolation_type") == 0)
        s.set_interpolation_type(static_cast<SIRFRegNiftyResample<float>::InterpolationType>(intDataFromHandle(hv)));
    else
        return parameterNotFound(name, __FILE__, __LINE__);
    return new DataHandle;
}
// get
void*
sirf::cSIRFReg_SIRFRegNiftyResampleParameter(const DataHandle* handle, const char* name)
{
    SIRFRegNiftyResample<float>& s = objectFromHandle<SIRFRegNiftyResample<float> >(handle);
    if (strcmp(name, "output") == 0)
        return newObjectHandle(s.get_output());
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}

// ------------------------------------------------------------------------------------ //
//   SIRFRegImageWeightedMean
// ------------------------------------------------------------------------------------ //
// get
void*
sirf::cSIRFReg_SIRFRegImageWeightedMeanParameter(const DataHandle* handle, const char* name)
{
    SIRFRegImageWeightedMean<float>& s = objectFromHandle<SIRFRegImageWeightedMean<float> >(handle);
    if (strcmp(name, "output") == 0)
        return newObjectHandle(s.get_output());
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}

// ------------------------------------------------------------------------------------ //
//   SIRFRegAffineTransformation
// ------------------------------------------------------------------------------------ //
// get
void*
sirf::cSIRFReg_SIRFRegAffineTransformationParameter(const DataHandle* handle, const char* name)
{
    SIRFRegAffineTransformation<float>& s = objectFromHandle<SIRFRegAffineTransformation<float> >(handle);
    if (strcmp(name, "determinant") == 0) {
        return dataHandle<float>(s.get_determinant());
    }
    if (strcmp(name,"identity") == 0) {
        std::shared_ptr<SIRFRegAffineTransformation<float> > sptr_id(new SIRFRegAffineTransformation<float>(SIRFRegAffineTransformation<float>::get_identity()));
        return newObjectHandle(sptr_id);
    }
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}
