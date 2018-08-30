/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2018 University College London.

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
#include "stir_x.h"
#include "csirfreg.h"
#include "csirfreg_p.h"
#include "SIRFRegMisc.h"
#include "SIRFImageData.h"
#include "SIRFImageDataDeformation.h"
#include "SIRFRegNiftyAladinSym.h"
#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFRegImageWeightedMean.h"
#include "stir_data_containers.h"

using namespace stir;
using namespace sirf;

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

//default constructors
extern "C"
void* cSIRFReg_newObject(const char* name)
{
	try {
		if (boost::iequals(name, "SIRFImageData"))
			return newObjectHandle<SIRFImageData>();
		if (boost::iequals(name, "SIRFImageDataDeformation"))
			return newObjectHandle<SIRFImageDataDeformation>();
        if (boost::iequals(name, "SIRFRegNiftyAladinSym"))
            return newObjectHandle<SIRFRegNiftyAladinSym<float> >();
        if (boost::iequals(name, "SIRFRegNiftyF3dSym"))
            return newObjectHandle<SIRFRegNiftyF3dSym<float> >();
        if (boost::iequals(name, "SIRFRegNiftyResample"))
            return newObjectHandle<SIRFRegNiftyResample>();
        if (boost::iequals(name, "SIRFRegImageWeightedMean3D"))
            return newObjectHandle<SIRFRegImageWeightedMean3D>();
        if (boost::iequals(name, "SIRFRegImageWeightedMean4D"))
            return newObjectHandle<SIRFRegImageWeightedMean4D>();
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

// set parameters
extern "C"
void* cSIRFReg_setParameter
(void* ptr_s, const char* obj, const char* name, const void* ptr_v)
{
	try {
		CAST_PTR(DataHandle, hs, ptr_s);
		CAST_PTR(DataHandle, hv, ptr_v);
        if (boost::iequals(obj, "SIRFReg"))
            return cSIRFReg_setSIRFRegParameter(ptr_s, name, ptr_v);
        if (boost::iequals(obj, "SIRFRegNiftyF3dSym"))
            return cSIRFReg_setSIRFRegNiftyF3dSymParameter(ptr_s, name, ptr_v);
        if (boost::iequals(obj, "SIRFRegNiftyResample"))
            return cSIRFReg_setSIRFRegNiftyResampleParameter(ptr_s, name, ptr_v);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

// get functions
extern "C"
void* cSIRFReg_parameter(const void* ptr, const char* obj, const char* name) 
{
	try {
		CAST_PTR(DataHandle, handle, ptr);
        if (boost::iequals(obj, "SIRFImageData"))
            return cSIRFReg_SIRFImageDataParameter(handle, name);
        if (boost::iequals(obj, "SIRFReg"))
            return cSIRFReg_SIRFRegParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegNiftyResample"))
            return cSIRFReg_SIRFRegNiftyResampleParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegImageWeightedMean3D"))
            return cSIRFReg_SIRFRegImageWeightedMean3DParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegImageWeightedMean4D"))
            return cSIRFReg_SIRFRegImageWeightedMean4DParameter(handle, name);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

// constructors from file
extern "C"
void* cSIRFReg_objectFromFile(const char* name, const char* filename)
{
	try {
		if (boost::iequals(name, "SIRFImageData")) {
			shared_ptr<SIRFImageData> 
				sptr(new SIRFImageData(filename));
			return newObjectHandle(sptr);
		}
		if (boost::iequals(name, "SIRFImageDataDeformation")) {
			shared_ptr<SIRFImageDataDeformation> 
				sptr(new SIRFImageDataDeformation(filename));
			return newObjectHandle(sptr);
		}
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

// -------------------------------------------------------------------------------- //
//      SIRFRegMisc
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_do_nifti_images_match(const void* im1, const void* im2, const float accuracy_percentage_of_max)
{
    try {
        SIRFImageData& img1 = objectFromHandle<SIRFImageData>(im1);
        SIRFImageData& img2 = objectFromHandle<SIRFImageData>(im2);
        return dataHandle(SIRFRegMisc::do_nifti_images_match(img1,img2,accuracy_percentage_of_max));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_filename(const char* filename)
{
    try {
        SIRFRegMisc::dump_nifti_info(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im1(const void* im1)
{
    try {
        SIRFImageData& img1 = objectFromHandle<SIRFImageData>(im1);
        SIRFRegMisc::dump_nifti_info(img1);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im2(const void* im1, const void* im2)
{
    try {
        std::vector<SIRFImageData> vec;
        vec.push_back(objectFromHandle<SIRFImageData>(im1));
        vec.push_back(objectFromHandle<SIRFImageData>(im2));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im3(const void* im1, const void* im2, const void* im3)
{
    try {
        std::vector<SIRFImageData> vec;
        vec.push_back(objectFromHandle<SIRFImageData>(im1));
        vec.push_back(objectFromHandle<SIRFImageData>(im2));
        vec.push_back(objectFromHandle<SIRFImageData>(im3));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im4(const void* im1, const void* im2, const void* im3, const void* im4)
{
    try {
        std::vector<SIRFImageData> vec;
        vec.push_back(objectFromHandle<SIRFImageData>(im1));
        vec.push_back(objectFromHandle<SIRFImageData>(im2));
        vec.push_back(objectFromHandle<SIRFImageData>(im3));
        vec.push_back(objectFromHandle<SIRFImageData>(im4));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im5(const void* im1, const void* im2, const void* im3, const void* im4, const void* im5)
{
    try {
        std::vector<SIRFImageData> vec;
        vec.push_back(objectFromHandle<SIRFImageData>(im1));
        vec.push_back(objectFromHandle<SIRFImageData>(im2));
        vec.push_back(objectFromHandle<SIRFImageData>(im3));
        vec.push_back(objectFromHandle<SIRFImageData>(im4));
        vec.push_back(objectFromHandle<SIRFImageData>(im5));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFImageData
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFImageData_from_PETImageData(void* ptr)
{
	try {
		sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(ptr);
		shared_ptr<SIRFImageData> 
			sptr(new SIRFImageData(pet_im));
		return newObjectHandle(sptr);
	}
	CATCH;
}
extern "C"
void* cSIRFReg_SIRFImageData_save_to_file(const void* ptr, const char* filename)
{
    try {
        SIRFImageData& im = objectFromHandle<SIRFImageData>(ptr);
        im.save_to_file(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFImageData_copy_data_to(const void* ptr, const void* obj)
{
    try {
        SIRFImageData& im = objectFromHandle<SIRFImageData>(ptr);
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(obj);
        im.copy_data_to(pet_im);
        return new DataHandle;
    }
    CATCH;
}

extern "C"
void* cSIRFReg_SIRFImageData_fill(const void* ptr, const float val)
{
    try {
        SIRFImageData& im = objectFromHandle<SIRFImageData>(ptr);
        im.fill(val);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFImageDataDeformation
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFImageDataDeformation_save_to_file(const void *ptr, const char* filename, const bool split_xyz)
{
	try {
		SIRFImageDataDeformation& im = objectFromHandle<SIRFImageDataDeformation>(ptr);
		im.save_to_file(filename,split_xyz);
		return new DataHandle;
	}
	CATCH;
}
extern "C"
void* cSIRFReg_SIRFImageDataDeformation_create_from_3D_image(void* ptr, const void* obj)
{
	try {
		SIRFImageDataDeformation& im = objectFromHandle<SIRFImageDataDeformation>(ptr);
		SIRFImageData& im2 = objectFromHandle<SIRFImageData>(obj);
		im.create_from_3D_image(im2);
		return new DataHandle;
	}
	CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFReg
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFReg_update(void* ptr)
{
    try {
        SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
        reg.update();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_save_image(const void *ptr, const char* filename)
{
	try {
		SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
		reg.save_warped_image(filename);
		return new DataHandle;
	}
	CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_save_deformation_displacement_image(const void *ptr, const char* filename, const char* type, const bool split_xyz)
{
    try {
        SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
        if (strcmp(type,"fwrd_deformation"))
            reg.save_deformation_field_fwrd_image(filename, split_xyz);
        else if (strcmp(type,"back_deformation"))
            reg.save_deformation_field_back_image(filename, split_xyz);
        else if (strcmp(type,"fwrd_displacement"))
            reg.save_displacement_field_fwrd_image(filename, split_xyz);
        else if (strcmp(type,"back_displacement"))
            reg.save_displacement_field_back_image(filename, split_xyz);
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      SIRFRegNiftyAladinSym
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(const void *ptr, const char* filename, const char* dir)
{
    try {
        SIRFRegNiftyAladinSym<float>& reg = objectFromHandle<SIRFRegNiftyAladinSym<float> >(ptr);
        if (strcmp(dir,"fwrd"))
            reg.save_transformation_matrix_fwrd(filename);
        else if (strcmp(dir,"back"))
            reg.save_transformation_matrix_fwrd(filename);
        else
            throw std::runtime_error("only accept fwrd or back as argument to dir for saving transformation matrix");
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      SIRFRegNiftyResample
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegNiftyResample_update(void* ptr)
{
    try {
        SIRFRegNiftyResample& res = objectFromHandle<SIRFRegNiftyResample>(ptr);
        res.update();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegNiftyResample_save_resampled_image(void* ptr, const char* filename)
{
    try {
        SIRFRegNiftyResample& res = objectFromHandle<SIRFRegNiftyResample>(ptr);
        res.save_resampled_image(filename);
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      SIRFRegImageWeightedMean3D
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean3D_add_image(void* ptr, const void *obj, const float weight)
{
    try {
        SIRFRegImageWeightedMean3D& im_weight = objectFromHandle<SIRFRegImageWeightedMean3D>(ptr);
        SIRFImageData& im = objectFromHandle<SIRFImageData>(obj);
        im_weight.add_image(im,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean3D_add_image_filename(void* ptr, const char* filename, const float weight)
{
    try {
        SIRFRegImageWeightedMean3D& im_weight = objectFromHandle<SIRFRegImageWeightedMean3D>(ptr);
        im_weight.add_image(filename,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean3D_update(void* ptr)
{
    try {
        SIRFRegImageWeightedMean3D& im_weight = objectFromHandle<SIRFRegImageWeightedMean3D>(ptr);
        im_weight.update();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean3D_save_image_to_file(const void *ptr, const char* filename)
{
    try {
        SIRFRegImageWeightedMean3D& im_weight = objectFromHandle<SIRFRegImageWeightedMean3D>(ptr);
        im_weight.save_image_to_file(filename);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFRegImageWeightedMean3D
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean4D_add_image(void* ptr, const void *obj, const float weight)
{
    try {
        SIRFRegImageWeightedMean3D& im_weight = objectFromHandle<SIRFRegImageWeightedMean3D>(ptr);
        SIRFImageData& im = objectFromHandle<SIRFImageData>(obj);
        im_weight.add_image(im,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean4D_add_image_filename(void* ptr, const char* filename, const float weight)
{
    try {
        SIRFRegImageWeightedMean4D& im_weight = objectFromHandle<SIRFRegImageWeightedMean4D>(ptr);
        im_weight.add_image(filename,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean4D_update(void* ptr)
{
    try {
        SIRFRegImageWeightedMean4D& im_weight = objectFromHandle<SIRFRegImageWeightedMean4D>(ptr);
        im_weight.update();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean4D_save_image_to_file(const void *ptr, const char* filename)
{
    try {
        SIRFRegImageWeightedMean4D& im_weight = objectFromHandle<SIRFRegImageWeightedMean4D>(ptr);
        im_weight.save_image_to_file(filename);
        return new DataHandle;
    }
    CATCH;
}
