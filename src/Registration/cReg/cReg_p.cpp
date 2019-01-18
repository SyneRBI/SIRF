/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2017 - 2019 University College London

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
#include "sirf/cReg/cReg_p.h"
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/Registration.h"
#include "sirf/cReg/NiftyF3dSym.h"
#include "sirf/cReg/NiftyResample.h"
#include "sirf/cReg/ImageWeightedMean.h"
#include "sirf/cReg/AffineTransformation.h"

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
sirf::cReg_NiftiImageDataParameter(const DataHandle* handle, const char* name)
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
//   Registration
// ------------------------------------------------------------------------------------ //
// set
void*
sirf::cReg_setRegistrationParameter(void* hp, const char* name, const void* hv)
{
    std::shared_ptr<const ImageData> im_sptr;

    Registration<float>& s = objectFromHandle<Registration<float> >(hp);
    if (strcmp(name, "parameter_file") == 0)
        s.set_parameter_file(charDataFromHandle(hv));
    else if (strcmp(name, "reference_image") == 0) {
        getObjectSptrFromHandle<const ImageData>(hv, im_sptr);
        s.set_reference_image(im_sptr);
    }
    else if (strcmp(name, "floating_image") == 0) {
        getObjectSptrFromHandle<const ImageData>(hv, im_sptr);
        std::cout << "\n setting floating image.\n";
        s.set_floating_image(im_sptr);
    }
    else if (strcmp(name, "reference_mask") == 0) {
        getObjectSptrFromHandle<const ImageData>(hv, im_sptr);
        s.set_reference_mask(im_sptr);
    }
    else if (strcmp(name, "floating_mask") == 0) {
        getObjectSptrFromHandle<const ImageData>(hv, im_sptr);
        s.set_floating_mask(im_sptr);
    }
	else
		return parameterNotFound(name, __FILE__, __LINE__);
	return new DataHandle;
}
// get
void*
sirf::cReg_RegistrationParameter(const DataHandle* handle, const char* name)
{
    Registration<float>& s = objectFromHandle<Registration<float> >(handle);
    if (strcmp(name, "output") == 0) {
        return newObjectHandle(s.get_output_sptr());
	}
	else
		return parameterNotFound(name, __FILE__, __LINE__);
}
// ------------------------------------------------------------------------------------ //
//   NiftyF3dSym
// ------------------------------------------------------------------------------------ //
// set
void*
sirf::cReg_setNiftyF3dSymParameter(void* hp, const char* name, const void* hv)
{
    NiftyF3dSym<float>& s = objectFromHandle<NiftyF3dSym<float> >(hp);
    if (strcmp(name, "floating_time_point") == 0)
        s.set_floating_time_point(intDataFromHandle(hv));
    else if (strcmp(name, "reference_time_point") == 0)
        s.set_reference_time_point(intDataFromHandle(hv));
    else if (strcmp(name, "initial_affine_transformation") == 0) {
        const std::shared_ptr<const AffineTransformation<float> > mat_sptr =
                std::make_shared<const AffineTransformation<float> >(objectFromHandle<const AffineTransformation<float> >(hv));
        s.set_initial_affine_transformation(mat_sptr);
    }
    else
        return parameterNotFound(name, __FILE__, __LINE__);
    return new DataHandle;
}
// ------------------------------------------------------------------------------------ //
//   NiftyResample
// ------------------------------------------------------------------------------------ //
// set
void*
sirf::cReg_setNiftyResampleParameter(void* hp, const char* name, const void* hv)
{
    std::shared_ptr<const ImageData> im_sptr;

    NiftyResample<float>& s = objectFromHandle<NiftyResample<float> >(hp);
    if (strcmp(name, "reference_image") == 0) {
        getObjectSptrFromHandle<const ImageData>(hv, im_sptr);
        s.set_reference_image(im_sptr);
    }
    else if (strcmp(name, "floating_image") == 0) {
        getObjectSptrFromHandle<const ImageData>(hv, im_sptr);
        s.set_floating_image(im_sptr);
    }
    else if (strcmp(name, "interpolation_type") == 0)
        s.set_interpolation_type(static_cast<NiftyResample<float>::InterpolationType>(intDataFromHandle(hv)));
    else
        return parameterNotFound(name, __FILE__, __LINE__);
    return new DataHandle;
}
// get
void*
sirf::cReg_NiftyResampleParameter(const DataHandle* handle, const char* name)
{
    NiftyResample<float>& s = objectFromHandle<NiftyResample<float> >(handle);
    if (strcmp(name, "output") == 0)
        return newObjectHandle(s.get_output_sptr());
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}

// ------------------------------------------------------------------------------------ //
//   ImageWeightedMean
// ------------------------------------------------------------------------------------ //
// get
void*
sirf::cReg_ImageWeightedMeanParameter(const DataHandle* handle, const char* name)
{
    ImageWeightedMean<float>& s = objectFromHandle<ImageWeightedMean<float> >(handle);
    if (strcmp(name, "output") == 0)
        return newObjectHandle(s.get_output_sptr());
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}

// ------------------------------------------------------------------------------------ //
//   AffineTransformation
// ------------------------------------------------------------------------------------ //
// get
void*
sirf::cReg_AffineTransformationParameter(const DataHandle* handle, const char* name)
{
    AffineTransformation<float>& s = objectFromHandle<AffineTransformation<float> >(handle);
    if (strcmp(name, "determinant") == 0) {
        return dataHandle<float>(s.get_determinant());
    }
    if (strcmp(name,"identity") == 0) {
        std::shared_ptr<AffineTransformation<float> > sptr_id(new AffineTransformation<float>);
        return newObjectHandle(sptr_id);
    }
    else
        return parameterNotFound(name, __FILE__, __LINE__);
}
