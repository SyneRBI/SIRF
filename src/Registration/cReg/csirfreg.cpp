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
#include "SIRFRegTransformation.h"
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
        if (boost::iequals(name, "SIRFRegTransformationAffine"))
            return newObjectHandle<SIRFRegTransformationAffine>();
        if (boost::iequals(name, "SIRFRegTransformationDisplacement"))
            return newObjectHandle<SIRFRegTransformationDisplacement>();
        if (boost::iequals(name, "SIRFRegTransformationDeformation"))
            return newObjectHandle<SIRFRegTransformationDeformation>();
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
        if (boost::iequals(name, "SIRFRegTransformationAffine")) {
            shared_ptr<SIRFRegTransformationAffine>
                sptr(new SIRFRegTransformationAffine(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "SIRFRegTransformationDisplacement")) {
            shared_ptr<SIRFRegTransformationDisplacement>
                sptr(new SIRFRegTransformationDisplacement(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "SIRFRegTransformationDeformation")) {
            shared_ptr<SIRFRegTransformationDeformation>
                sptr(new SIRFRegTransformationDeformation(filename));
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
extern "C"
void* cSIRFReg_SIRFReg_open_TM(const char* filename, size_t ptr_TM)
{
    try {
        float* TM = (float*)ptr_TM;
        mat44 tm;
        SIRFRegMisc::open_transformation_matrix(tm, filename);

        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                TM[i+j*4] = tm.m[i][j];
        return new DataHandle;
    }
    CATCH;
}

void* cSIRFReg_compose_transformations_into_single_deformation2(const void* im, const void* trans1, const void* trans2)
{
    try {
        SIRFImageData& ref = objectFromHandle<SIRFImageData>(im);
        shared_ptr<SIRFRegTransformationDeformation> def_sptr(new SIRFRegTransformationDeformation());
        std::vector<SIRFRegTransformation*> vec;
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans1));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans2));

        SIRFRegMisc::compose_transformations_into_single_deformation(*def_sptr, vec, ref);
        return newObjectHandle(def_sptr);
    }
    CATCH;
}
void* cSIRFReg_compose_transformations_into_single_deformation3(const void* im, const void* trans1, const void* trans2, const void* trans3)
{
    try {
        SIRFImageData& ref = objectFromHandle<SIRFImageData>(im);
        shared_ptr<SIRFRegTransformationDeformation> def_sptr(new SIRFRegTransformationDeformation());
        std::vector<SIRFRegTransformation*> vec;
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans1));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans2));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans3));

        SIRFRegMisc::compose_transformations_into_single_deformation(*def_sptr, vec, ref);
        return newObjectHandle(def_sptr);
    }
    CATCH;
}
void* cSIRFReg_compose_transformations_into_single_deformation4(const void* im, const void* trans1, const void* trans2, const void* trans3, const void* trans4)
{
    try {
        SIRFImageData& ref = objectFromHandle<SIRFImageData>(im);
        shared_ptr<SIRFRegTransformationDeformation> def_sptr(new SIRFRegTransformationDeformation());
        std::vector<SIRFRegTransformation*> vec;
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans1));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans2));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans3));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans4));

        SIRFRegMisc::compose_transformations_into_single_deformation(*def_sptr, vec, ref);
        return newObjectHandle(def_sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_compose_transformations_into_single_deformation5(const void* im, const void* trans1, const void* trans2, const void* trans3, const void* trans4, const void* trans5)
{
    try {
        SIRFImageData& ref = objectFromHandle<SIRFImageData>(im);
        shared_ptr<SIRFRegTransformationDeformation> def_sptr(new SIRFRegTransformationDeformation());
        std::vector<SIRFRegTransformation*> vec;
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans1));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans2));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans3));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans4));
        vec.push_back(&objectFromHandle<SIRFRegTransformation>(trans5));

        SIRFRegMisc::compose_transformations_into_single_deformation(*def_sptr, vec, ref);
        return newObjectHandle(def_sptr);
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
extern "C"
void* cSIRFReg_SIRFImageData_deep_copy(const void* ptr)
{
    try {
        SIRFImageData& im = objectFromHandle<SIRFImageData>(ptr);
        shared_ptr<SIRFImageData> sptr(new SIRFImageData(im.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFImageData_get_dimensions(const void* ptr, size_t ptr_dim)
{
    try {
        SIRFImageData& im = objectFromHandle<SIRFImageData>(ptr);
        int* dim = (int*)ptr_dim;
        im.get_dimensions(dim);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFImageData_get_data(const void* ptr, size_t ptr_data)
{
    try {
        SIRFImageData& im = objectFromHandle<SIRFImageData>(ptr);
        SIRFImageData copy = im.deep_copy();
        SIRFRegMisc::change_datatype<float>(copy);
        float* data = (float*)ptr_data;
        size_t mem = im.get_raw_nifti_sptr()->nvox * size_t(im.get_raw_nifti_sptr()->nbyper);
        // Copy!
        memcpy(data, im.get_raw_nifti_sptr()->data, mem);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFImageData_maths(const void* ptr, const void* obj, const int maths_type)
{
    try {
        SIRFImageData& im1 = objectFromHandle<SIRFImageData>(ptr);
        SIRFImageData& im2 = objectFromHandle<SIRFImageData>(obj);
        SIRFImageData res;
        if (maths_type == 1)
            res = im1 + im2;
        else if (maths_type == -1)
            res = im1 - im2;
        shared_ptr<SIRFImageData> sptr(new SIRFImageData(res));
        return newObjectHandle(sptr);
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
extern "C"
void* cSIRFReg_SIRFImageDataDeformation_deep_copy(const void* ptr)
{
    try {
        SIRFImageDataDeformation& im = objectFromHandle<SIRFImageDataDeformation>(ptr);
        shared_ptr<SIRFImageDataDeformation> sptr(new SIRFImageDataDeformation(im.deep_copy()));
        return newObjectHandle(sptr);
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
void* cSIRFReg_SIRFReg_get_deformation_displacement_image(const void* ptr, const char* type)
{
    try {
        SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
        shared_ptr<SIRFImageDataDeformation> sptr;
        if (strcmp(type,"fwrd_deformation") == 0)
            sptr.reset(new SIRFImageDataDeformation(reg.get_deformation_field_fwrd()));
        else if (strcmp(type,"back_deformation") == 0)
            sptr.reset(new SIRFImageDataDeformation(reg.get_deformation_field_back()));
        else if (strcmp(type,"fwrd_displacement") == 0)
            sptr.reset(new SIRFImageDataDeformation(reg.get_displacement_field_fwrd()));
        else if (strcmp(type,"back_displacement") == 0)
            sptr.reset(new SIRFImageDataDeformation(reg.get_displacement_field_back()));
        else
            throw std::runtime_error("cSIRFReg_SIRFReg_get_deformation_displacement_image: Bad return type.");
        return newObjectHandle(sptr);
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
        if (strcmp(dir,"fwrd") == 0)
            reg.save_transformation_matrix_fwrd(filename);
        else if (strcmp(dir,"back") == 0)
            reg.save_transformation_matrix_fwrd(filename);
        else
            throw std::runtime_error("only accept fwrd or back as argument to dir for saving transformation matrix");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_get_TM(const void* ptr, size_t ptr_TM, const char* dir)
{
    try {
        SIRFRegNiftyAladinSym<float>& reg = objectFromHandle<SIRFRegNiftyAladinSym<float> >(ptr);
        float* TM = (float*)ptr_TM;
        mat44 trans_m;
        if (strcmp(dir,"fwrd") == 0)
            trans_m = reg.get_transformation_matrix_fwrd();
        else if (strcmp(dir,"back") == 0)
            trans_m = reg.get_transformation_matrix_back();
        else
            throw std::runtime_error("only accept fwrd or back as argument to dir for saving transformation matrix");

        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                TM[i+j*4] = trans_m.m[i][j];

        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      SIRFRegNiftyResample
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegNiftyResample_add_transformation(void* self, const void* trans, const char *type)
{
    try {
        SIRFRegNiftyResample& res = objectFromHandle<SIRFRegNiftyResample>(self);
        if (strcmp(type, "affine") == 0)
            res.add_transformation_affine(objectFromHandle<SIRFRegTransformationAffine>(trans));
        else if (strcmp(type, "displacement") == 0)
            res.add_transformation_disp(objectFromHandle<SIRFRegTransformationDisplacement>(trans));
        else if (strcmp(type, "deformation") == 0)
            res.add_transformation_def(objectFromHandle<SIRFRegTransformationDeformation>(trans));
        else
            throw std::runtime_error("only accept 'affine', 'displacement' or 'deformation' as argument adding transformation matrix to resample.");
        return new DataHandle;
    }
    CATCH;
}
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

// -------------------------------------------------------------------------------- //
//      SIRFRegTransformation
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegTransformation_get_as_deformation_field(const void* ptr, const void* ref)
{
    try {
        SIRFRegTransformation& trans = objectFromHandle<SIRFRegTransformation>(ptr);
        SIRFImageData& ref_im = objectFromHandle<SIRFImageData>(ref);
        SIRFImageDataDeformation res = trans.get_as_deformation_field(ref_im);
        shared_ptr<SIRFImageData> sptr(new SIRFImageData(res));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegTransformationAffine_construct_from_TM(size_t ptr_TM)
{
    try {
        float* TM = (float*)ptr_TM;
        mat44 trans_m;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                trans_m.m[i][j] = TM[i+j*4];

        SIRFRegTransformationAffine trans(trans_m);
        shared_ptr<SIRFRegTransformationAffine> sptr(new SIRFRegTransformationAffine(trans));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegTransformationDisplacement_construct_from_SIRFImageDataDeformation(const void* ptr)
{
    try {
        SIRFImageDataDeformation& disp = objectFromHandle<SIRFImageDataDeformation>(ptr);
        SIRFRegTransformationDisplacement trans(disp);
        shared_ptr<SIRFRegTransformationDisplacement> sptr(new SIRFRegTransformationDisplacement(trans.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegTransformationDeformation_construct_from_SIRFImageDataDeformation(const void* ptr)
{
    try {
        SIRFImageDataDeformation& def = objectFromHandle<SIRFImageDataDeformation>(ptr);
        SIRFRegTransformationDeformation trans(def);
        shared_ptr<SIRFRegTransformationDeformation> sptr(new SIRFRegTransformationDeformation(trans.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
