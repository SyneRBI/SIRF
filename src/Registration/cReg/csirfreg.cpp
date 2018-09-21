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
#include "NiftiImage3D.h"
#include "NiftiImage3DTensor.h"
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
        if (boost::iequals(name, "NiftiImage"))
            return newObjectHandle<NiftiImage>();
        if (boost::iequals(name, "NiftiImage3D"))
            return newObjectHandle<NiftiImage3D>();
        if (boost::iequals(name, "NiftiImage3DTensor"))
            return newObjectHandle<NiftiImage3DTensor>();
        if (boost::iequals(name, "NiftiImage3DDisplacement"))
            return newObjectHandle<NiftiImage3DDisplacement>();
        if (boost::iequals(name, "NiftiImage3DDeformation"))
            return newObjectHandle<NiftiImage3DDeformation>();
        if (boost::iequals(name, "SIRFRegNiftyAladinSym"))
            return newObjectHandle<SIRFRegNiftyAladinSym<float> >();
        if (boost::iequals(name, "SIRFRegNiftyF3dSym"))
            return newObjectHandle<SIRFRegNiftyF3dSym<float> >();
        if (boost::iequals(name, "SIRFRegNiftyResample"))
            return newObjectHandle<SIRFRegNiftyResample>();
        if (boost::iequals(name, "SIRFRegImageWeightedMean"))
            return newObjectHandle<SIRFRegImageWeightedMean>();
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
        if (boost::iequals(obj, "NiftiImage"))
            return cSIRFReg_NiftiImageParameter(handle, name);
        if (boost::iequals(obj, "SIRFReg"))
            return cSIRFReg_SIRFRegParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegNiftyResample"))
            return cSIRFReg_SIRFRegNiftyResampleParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegImageWeightedMean"))
            return cSIRFReg_SIRFRegImageWeightedMeanParameter(handle, name);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

// constructors from file
extern "C"
void* cSIRFReg_objectFromFile(const char* name, const char* filename)
{
	try {
        if (boost::iequals(name, "NiftiImage")) {
            shared_ptr<NiftiImage>
                sptr(new NiftiImage(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImage3D")) {
            shared_ptr<NiftiImage3D>
                sptr(new NiftiImage3D(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImage3DTensor")) {
            shared_ptr<NiftiImage3DTensor>
                sptr(new NiftiImage3DTensor(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImage3DDisplacement")) {
            shared_ptr<NiftiImage3DDisplacement>
                sptr(new NiftiImage3DDisplacement(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImage3DDeformation")) {
            shared_ptr<NiftiImage3DDeformation>
                sptr(new NiftiImage3DDeformation(filename));
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
        NiftiImage& img1 = objectFromHandle<NiftiImage>(im1);
        NiftiImage& img2 = objectFromHandle<NiftiImage>(im2);
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
        NiftiImage& img1 = objectFromHandle<NiftiImage>(im1);
        SIRFRegMisc::dump_nifti_info(img1);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im2(const void* im1, const void* im2)
{
    try {
        std::vector<NiftiImage> vec;
        vec.push_back(objectFromHandle<NiftiImage>(im1));
        vec.push_back(objectFromHandle<NiftiImage>(im2));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im3(const void* im1, const void* im2, const void* im3)
{
    try {
        std::vector<NiftiImage> vec;
        vec.push_back(objectFromHandle<NiftiImage>(im1));
        vec.push_back(objectFromHandle<NiftiImage>(im2));
        vec.push_back(objectFromHandle<NiftiImage>(im3));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im4(const void* im1, const void* im2, const void* im3, const void* im4)
{
    try {
        std::vector<NiftiImage> vec;
        vec.push_back(objectFromHandle<NiftiImage>(im1));
        vec.push_back(objectFromHandle<NiftiImage>(im2));
        vec.push_back(objectFromHandle<NiftiImage>(im3));
        vec.push_back(objectFromHandle<NiftiImage>(im4));
        SIRFRegMisc::dump_nifti_info(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_dump_nifti_info_im5(const void* im1, const void* im2, const void* im3, const void* im4, const void* im5)
{
    try {
        std::vector<NiftiImage> vec;
        vec.push_back(objectFromHandle<NiftiImage>(im1));
        vec.push_back(objectFromHandle<NiftiImage>(im2));
        vec.push_back(objectFromHandle<NiftiImage>(im3));
        vec.push_back(objectFromHandle<NiftiImage>(im4));
        vec.push_back(objectFromHandle<NiftiImage>(im5));
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
        NiftiImage3D& ref = objectFromHandle<NiftiImage3D>(im);
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
        NiftiImage3D& ref = objectFromHandle<NiftiImage3D>(im);
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
        NiftiImage3D& ref = objectFromHandle<NiftiImage3D>(im);
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
        NiftiImage3D& ref = objectFromHandle<NiftiImage3D>(im);
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
//      NiftiImage
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImage_save_to_file(const void* ptr, const char* filename)
{
	try {
        NiftiImage& im = objectFromHandle<NiftiImage>(ptr);
        im.save_to_file(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_fill(const void* ptr, const float val)
{
    try {
        NiftiImage& im = objectFromHandle<NiftiImage>(ptr);
        im.fill(val);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_deep_copy(const void* copy_ptr, const void* orig_ptr)
{
    try {
        NiftiImage& orig = objectFromHandle<NiftiImage>(orig_ptr);
        NiftiImage& copy = objectFromHandle<NiftiImage>(copy_ptr);
        copy = orig.deep_copy();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_get_dimensions(const void* ptr, size_t ptr_dim)
{
    try {
        NiftiImage& im = objectFromHandle<NiftiImage>(ptr);
        int* dim = (int*)ptr_dim;
        im.get_dimensions(dim);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_get_data(const void* ptr, size_t ptr_data)
{
    try {
        NiftiImage& im = objectFromHandle<NiftiImage>(ptr);
        NiftiImage copy = im.deep_copy();
        SIRFRegMisc::change_datatype<float>(copy);
        float* data = (float*)ptr_data;
        size_t mem = copy.get_raw_nifti_sptr()->nvox * size_t(copy.get_raw_nifti_sptr()->nbyper);
        // Copy!
        memcpy(data, copy.get_raw_nifti_sptr()->data, mem);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_maths(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type)
{
    try {
        if (abs(maths_type) != 1)
            throw std::runtime_error("cSIRFReg_NiftiImage_maths: Bad maths type (1: addition, -1: subtraction.");

        NiftiImage& res = objectFromHandle<NiftiImage>(res_ptr);
        NiftiImage& im1 = objectFromHandle<NiftiImage>(im1_ptr);
        NiftiImage& im2 = objectFromHandle<NiftiImage>(im2_ptr);

        if (maths_type ==  1) res = im1 + im2;
        if (maths_type == -1) res = im1 - im2;
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImage3D
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImage3D_from_PETImageData(void* ptr)
{
	try {
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(ptr);
        shared_ptr<NiftiImage3D>
            sptr(new NiftiImage3D(pet_im));
        return newObjectHandle(sptr);
    }
	CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage3D_copy_data_to(const void* ptr, const void* obj)
{
    try {
        NiftiImage3D& im = objectFromHandle<NiftiImage3D>(ptr);
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(obj);
        im.copy_data_to(pet_im);
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      NiftiImage3DTensor
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImage3DTensor_save_to_file_split_xyz_components(const void *ptr, const char* filename)
{
	try {
        NiftiImage3DTensor& im = objectFromHandle<NiftiImage3DTensor>(ptr);
        im.save_to_file_split_xyz_components(filename);
		return new DataHandle;
	}
	CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage3DTensor_create_from_3D_image(const void* ptr, const void* obj)
{
    try {
        NiftiImage3DTensor& im = objectFromHandle<NiftiImage3DTensor>(ptr);
        NiftiImage3D& im3d = objectFromHandle<NiftiImage3D>(obj);
        im.create_from_3D_image(im3d);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr)
{
    try {
        NiftiImage3D& x = objectFromHandle<NiftiImage3D>(x_ptr);
        NiftiImage3D& y = objectFromHandle<NiftiImage3D>(y_ptr);
        NiftiImage3D& z = objectFromHandle<NiftiImage3D>(z_ptr);

        shared_ptr<NiftiImage3DTensor> sptr;
        if (strcmp(obj,"NiftiImage3DTensor"))
            sptr.reset(new NiftiImage3DTensor(x,y,z));
        else if (strcmp(obj,"NiftiImage3DDisplacement"))
            sptr.reset(new NiftiImage3DDisplacement(x,y,z));
        else if (strcmp(obj,"NiftiImage3DDeformation"))
            sptr.reset(new NiftiImage3DDeformation(x,y,z));
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
void* cSIRFReg_SIRFReg_get_deformation_displacement_image(const void* ptr, const char *transform_type)
{
    try {
        SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
        shared_ptr<NiftiImage3DDeformation> sptr;
        if (strcmp(transform_type, "fwrd_deformation") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDeformation>(new NiftiImage3DDeformation(reg.get_deformation_field_fwrd())));
        else if (strcmp(transform_type, "back_deformation") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDeformation>(new NiftiImage3DDeformation(reg.get_deformation_field_back())));
        else if (strcmp(transform_type, "fwrd_displacement") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDisplacement>(new NiftiImage3DDisplacement(reg.get_displacement_field_fwrd())));
        else if (strcmp(transform_type, "back_displacement") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDisplacement>(new NiftiImage3DDisplacement(reg.get_displacement_field_back())));
        else
            throw std::runtime_error("cSIRFReg_SIRFReg_get_deformation_displacement_image: Bad return type.");
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFRegNiftyAladinSym
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegNiftyAladinSym_save_transformation_matrix(const void *ptr, const char* filename, const char *dir)
{
    try {
        SIRFRegNiftyAladinSym<float>& reg = objectFromHandle<SIRFRegNiftyAladinSym<float> >(ptr);
        if (strcmp(dir, "fwrd") == 0)
            reg.save_transformation_matrix_fwrd(filename);
        else if (strcmp(dir, "back") == 0)
            reg.save_transformation_matrix_fwrd(filename);
        else
            throw std::runtime_error("only accept fwrd or back as argument to dir for saving transformation matrix");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_get_TM(const void* ptr, size_t ptr_TM, const char *dir)
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
//      SIRFRegImageWeightedMean
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean_add_image(void* ptr, const void *obj, const float weight)
{
    try {
        SIRFRegImageWeightedMean& im_weight = objectFromHandle<SIRFRegImageWeightedMean>(ptr);
        NiftiImage& im = objectFromHandle<NiftiImage>(obj);
        im_weight.add_image(im,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight)
{
    try {
        SIRFRegImageWeightedMean& im_weight = objectFromHandle<SIRFRegImageWeightedMean>(ptr);
        im_weight.add_image(filename,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean_update(void* ptr)
{
    try {
        SIRFRegImageWeightedMean& im_weight = objectFromHandle<SIRFRegImageWeightedMean>(ptr);
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
        NiftiImage& ref_im = objectFromHandle<NiftiImage>(ref);
        NiftiImage3DDeformation res = trans.get_as_deformation_field(ref_im);
        shared_ptr<NiftiImage> sptr(new NiftiImage(res));
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
void* cSIRFReg_SIRFRegTransformationDisplacement_construct_from_NiftiImage3DDisplacement(const void* ptr)
{
    try {
        NiftiImage3DDisplacement& disp = objectFromHandle<NiftiImage3DDisplacement>(ptr);
        SIRFRegTransformationDisplacement trans(disp);
        shared_ptr<SIRFRegTransformationDisplacement> sptr(new SIRFRegTransformationDisplacement(trans.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegTransformationDeformation_construct_from_NiftiImage3DDeformation(const void* ptr)
{
    try {
        NiftiImage3DDeformation& def = objectFromHandle<NiftiImage3DDeformation>(ptr);
        SIRFRegTransformationDeformation trans(def);
        shared_ptr<SIRFRegTransformationDeformation> sptr(new SIRFRegTransformationDeformation(trans.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
