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
#include "SIRFRegMat44.h"
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
        if (boost::iequals(name, "SIRFRegMat44"))
            return newObjectHandle<SIRFRegMat44>();
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
        if (boost::iequals(obj, "SIRFRegMat44"))
            return cSIRFReg_SIRFRegMat44Parameter(handle, name);
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
        if (boost::iequals(name, "SIRFRegMat44")) {
            shared_ptr<SIRFRegMat44>
                sptr(new SIRFRegMat44(filename));
            return newObjectHandle(sptr);
        }
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImage
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImage_print_headers(const int num_ims, const void* im1, const void* im2, const void* im3, const void* im4, const void* im5)
{
    try {
        std::vector<NiftiImage> vec;
        if (num_ims >= 1) vec.push_back(objectFromHandle<NiftiImage>(im1));
        if (num_ims >= 2) vec.push_back(objectFromHandle<NiftiImage>(im2));
        if (num_ims >= 3) vec.push_back(objectFromHandle<NiftiImage>(im3));
        if (num_ims >= 4) vec.push_back(objectFromHandle<NiftiImage>(im4));
        if (num_ims >= 5) vec.push_back(objectFromHandle<NiftiImage>(im5));
        NiftiImage::print_headers(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_save_to_file(const void* ptr, const char* filename, const int datatype)
{
	try {
        NiftiImage& im = objectFromHandle<NiftiImage>(ptr);
        im.save_to_file(filename,datatype);
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
        for (int i=0; i<8; ++i)
            dim[i] = im.get_dimensions()[i];
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
        float* data = (float*)ptr_data;
        size_t mem = copy.get_raw_nifti_sptr()->nvox * size_t(copy.get_raw_nifti_sptr()->nbyper);
        // Copy!
        memcpy(data, copy.get_raw_nifti_sptr()->data, mem);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type)
{
    try {
        NiftiImage& res = objectFromHandle<NiftiImage>(res_ptr);
        NiftiImage& im1 = objectFromHandle<NiftiImage>(im1_ptr);
        NiftiImage& im2 = objectFromHandle<NiftiImage>(im2_ptr);

        if      (maths_type == 0) res = im1 + im2;
        else if (maths_type == 1) res = im1 - im2;
        else
            throw std::runtime_error("cSIRFReg_NiftiImage_maths_im: Bad maths type (0=add, 1=subtract).");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_maths_num(const void *res_ptr, const void* im1_ptr, const float val, const int maths_type)
{
    try {
        NiftiImage& res = objectFromHandle<NiftiImage>(res_ptr);
        NiftiImage& im1 = objectFromHandle<NiftiImage>(im1_ptr);

        if      (maths_type == 0) res = im1 + val;
        else if (maths_type == 1) res = im1 - val;
        else if (maths_type == 2) res = im1 * val;
        else
            throw std::runtime_error("cSIRFReg_NiftiImage_maths_val: Bad maths type (0=add, 1=subtract, 2=multiply.");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_equal(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImage& im1 = objectFromHandle<NiftiImage>(im1_ptr);
        NiftiImage& im2 = objectFromHandle<NiftiImage>(im2_ptr);
        return dataHandle<int>(im1 == im2);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_norm(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImage& im1 = objectFromHandle<NiftiImage>(im1_ptr);
        NiftiImage& im2 = objectFromHandle<NiftiImage>(im2_ptr);
        return dataHandle<float>(im1.get_norm(im2));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_get_original_datatype(const void* im_ptr)
{
    try {
        NiftiImage& im = objectFromHandle<NiftiImage>(im_ptr);
        return charDataHandleFromCharData(nifti_datatype_to_string(im.get_original_datatype()));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage_crop(const void* im_ptr, size_t min_index_ptr, size_t max_index_ptr)
{
    try {
        NiftiImage& im = objectFromHandle<NiftiImage>(im_ptr);
        int* min_index = (int*)min_index_ptr;
        int* max_index = (int*)max_index_ptr;
        im.crop(min_index,max_index);
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
void* cSIRFReg_NiftiImage3DTensor_save_to_file_split_xyz_components(const void *ptr, const char* filename, const int datatype)
{
	try {
        NiftiImage3DTensor& im = objectFromHandle<NiftiImage3DTensor>(ptr);
        im.save_to_file_split_xyz_components(filename, datatype);
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
extern "C"
void* cSIRFReg_NiftiImage3DTensor_flip_component(const void *ptr, const int dim)
{
    try {
        NiftiImage3DTensor& im = objectFromHandle<NiftiImage3DTensor>(ptr);
        im.flip_component(dim);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImage3DDeformation
// -------------------------------------------------------------------------------- //
void* cSIRFReg_NiftiImage3DDeformation_compose_single_deformation(const void* im, const int num_elements, const char* types, const void* trans1, const void* trans2, const void* trans3, const void* trans4, const void* trans5)
{
    try {
        // This is an ugly hack because I can't get virtual methods to work for multiple inherited (NiftiImage3DDeformation/NiftiImage3DDisplacement).
        // So, we also give a string which tells us what type they are, and we change the template type of objectFromHandle accordingly.

        // Also, we can't have default arguments in C, so if we only want to compose 3 transformations, set the 4th and 5th as 'None' in Python. In C,
        // we create a vector from all the non-'None' arguments and then convert them to their derived classes.

        // Sorry this is so ugly.

        // There's always going to be at least two transformations, so start by putting them in the vector
        std::vector<const void*> vec = {trans1, trans2};
        // Add in any extras, depending on the number of transformations
        if (num_elements >= 3) vec.push_back(trans3);
        if (num_elements >= 4) vec.push_back(trans4);
        if (num_elements >= 5) vec.push_back(trans5);

        // Vector for casting to the correct type
        std::vector<SIRFRegTransformation*> trans_vec;
        for (int i=0; i<num_elements; ++i)
            if      (types[i] == '1')
                trans_vec.push_back(&objectFromHandle<SIRFRegMat44>(vec.at(i)));
            else if (types[i] == '2')
                trans_vec.push_back(&objectFromHandle<NiftiImage3DDisplacement>(vec.at(i)));
            else if (types[i] == '3')
                trans_vec.push_back(&objectFromHandle<NiftiImage3DDeformation>(vec.at(i)));

        NiftiImage3D& ref = objectFromHandle<NiftiImage3D>(im);
        shared_ptr<NiftiImage3DDeformation> def_sptr
                (new NiftiImage3DDeformation(NiftiImage3DDeformation::compose_single_deformation(trans_vec, ref).deep_copy()));
        return newObjectHandle(def_sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImage3DDeformation_create_from_disp(const void* ptr, const void* disp_ptr)
{
    try {
        NiftiImage3DDeformation&  def  = objectFromHandle<NiftiImage3DDeformation>(ptr);
        NiftiImage3DDisplacement& disp = objectFromHandle<NiftiImage3DDisplacement>(disp_ptr);
        def.create_from_disp(disp);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImage3DDisplacement
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImage3DDisplacement_create_from_def(const void* ptr, const void* def_ptr)
{
    try {
        NiftiImage3DDisplacement& disp = objectFromHandle<NiftiImage3DDisplacement>(ptr);
        NiftiImage3DDeformation&  def  = objectFromHandle<NiftiImage3DDeformation>(def_ptr);
        disp.create_from_def(def);
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      SIRFReg
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFReg_process(void* ptr)
{
    try {
        SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
        reg.process();
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
        if (strcmp(transform_type, "forward_deformation") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDeformation>(new NiftiImage3DDeformation(reg.get_deformation_field_forward())));
        else if (strcmp(transform_type, "inverse_deformation") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDeformation>(new NiftiImage3DDeformation(reg.get_deformation_field_inverse())));
        else if (strcmp(transform_type, "forward_displacement") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDisplacement>(new NiftiImage3DDisplacement(reg.get_displacement_field_forward())));
        else if (strcmp(transform_type, "inverse_displacement") == 0)
            return newObjectHandle(shared_ptr<NiftiImage3DDisplacement>(new NiftiImage3DDisplacement(reg.get_displacement_field_inverse())));
        else
            throw std::runtime_error("cSIRFReg_SIRFReg_get_deformation_displacement_image: Bad return type.");
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_set_parameter(const void* ptr, const char* par, const char* arg1, const char* arg2)
{
    try {
        SIRFReg& reg = objectFromHandle<SIRFReg>(ptr);
        reg.set_parameter(par, arg1, arg2);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFRegNiftyAladinSym
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFReg_get_TM(const void* ptr, const char* dir)
{
    try {
        SIRFRegNiftyAladinSym<float>& reg = objectFromHandle<SIRFRegNiftyAladinSym<float> >(ptr);
        shared_ptr<SIRFRegMat44> sptr;
        if (strcmp(dir, "forward") == 0)
            sptr.reset(new SIRFRegMat44(reg.get_transformation_matrix_forward().deep_copy()));
        else if (strcmp(dir, "inverse") == 0)
            sptr.reset(new SIRFRegMat44(reg.get_transformation_matrix_inverse().deep_copy()));
        else
            throw std::runtime_error("only accept forward or inverse as argument to dir for saving transformation matrix");
        return newObjectHandle(sptr);
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
            res.add_transformation_affine(objectFromHandle<SIRFRegMat44>(trans));
        else if (strcmp(type, "displacement") == 0)
            res.add_transformation_disp(objectFromHandle<NiftiImage3DDisplacement>(trans));
        else if (strcmp(type, "deformation") == 0)
            res.add_transformation_def(objectFromHandle<NiftiImage3DDeformation>(trans));
        else
            throw std::runtime_error("only accept 'affine', 'displacement' or 'deformation' as argument adding transformation matrix to resample.");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegNiftyResample_process(void* ptr)
{
    try {
        SIRFRegNiftyResample& res = objectFromHandle<SIRFRegNiftyResample>(ptr);
        res.process();
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
        im_weight.add_image(NiftiImage(filename),weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean_process(void* ptr)
{
    try {
        SIRFRegImageWeightedMean& im_weight = objectFromHandle<SIRFRegImageWeightedMean>(ptr);
        im_weight.process();
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFRegTransformation
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegTransformation_get_as_deformation_field(const void* ptr, const char* name, const void* ref)
{
    try {
        SIRFRegTransformation *trans;

        if (strcmp(name,"SIRFRegMat44") == 0)
            trans = &objectFromHandle<SIRFRegMat44>(ptr);
        else if (strcmp(name,"NiftiImage3DDisplacement") == 0)
            trans = &objectFromHandle<NiftiImage3DDisplacement>(ptr);
        else if (strcmp(name,"NiftiImage3DDeformation") == 0)
            trans = &objectFromHandle<NiftiImage3DDeformation>(ptr);
        else
            throw std::runtime_error("cSIRFReg_SIRFRegTransformation_get_as_deformation_field: type should be affine, disp or def.");

        NiftiImage& ref_im = objectFromHandle<NiftiImage>(ref);
        shared_ptr<NiftiImage3DDeformation> sptr
                (new NiftiImage3DDeformation(trans->get_as_deformation_field(ref_im)));

        return newObjectHandle(sptr);
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFRegMat44
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegMat44_construct_from_TM(size_t ptr_TM)
{
    try {
        float* TM = (float*)ptr_TM;
        mat44 trans_m;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                trans_m.m[i][j] = TM[i+j*4];

        SIRFRegMat44 trans(trans_m);
        shared_ptr<SIRFRegMat44> sptr(new SIRFRegMat44(trans));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegMat44_deep_copy(const void* ptr)
{
    try {
        SIRFRegMat44& mat = objectFromHandle<SIRFRegMat44>(ptr);
        shared_ptr<SIRFRegMat44> sptr(new SIRFRegMat44(mat.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegMat44_save_to_file(const void* ptr, const char* filename)
{
    try {
        SIRFRegMat44& mat = objectFromHandle<SIRFRegMat44>(ptr);
        mat.save_to_file(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegMat44_fill(const void* ptr, const float val)
{
    try {
        SIRFRegMat44& mat = objectFromHandle<SIRFRegMat44>(ptr);
        mat.fill(val);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegMat44_as_array(const void* ptr, size_t ptr_TM)
{
    try {
        SIRFRegMat44& tm = objectFromHandle<SIRFRegMat44>(ptr);
        float* TM = (float*)ptr_TM;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                TM[i*4+j] = tm[i][j];
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegMat44_get_identity()
{
    try {
        shared_ptr<SIRFRegMat44> sptr(new SIRFRegMat44(SIRFRegMat44::get_identity()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegMat44_mul(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        SIRFRegMat44& mat1 = objectFromHandle<SIRFRegMat44>(mat1_ptr);
        SIRFRegMat44& mat2 = objectFromHandle<SIRFRegMat44>(mat2_ptr);
        shared_ptr<SIRFRegMat44> sptr(new SIRFRegMat44(mat1*mat2));
        return newObjectHandle(sptr);
    }
    CATCH;
}

extern "C"
void* cSIRFReg_SIRFRegMat44_equal(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        SIRFRegMat44& mat1 = objectFromHandle<SIRFRegMat44>(mat1_ptr);
        SIRFRegMat44& mat2 = objectFromHandle<SIRFRegMat44>(mat2_ptr);
        return dataHandle<int>(mat1 == mat2);
    }
    CATCH;
}
