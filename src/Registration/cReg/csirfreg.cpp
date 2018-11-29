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
#include "NiftiImageData3D.h"
#include "NiftiImageData3DTensor.h"
#include "SIRFRegNiftyAladinSym.h"
#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFRegImageWeightedMean.h"
#include "SIRFRegTransformation.h"
#include "SIRFRegAffineTransformation.h"
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
        if (boost::iequals(name, "NiftiImageData"))
            return newObjectHandle<NiftiImageData>();
        if (boost::iequals(name, "NiftiImageData3D"))
            return newObjectHandle<NiftiImageData3D>();
        if (boost::iequals(name, "NiftiImageData3DTensor"))
            return newObjectHandle<NiftiImageData3DTensor>();
        if (boost::iequals(name, "NiftiImageData3DDisplacement"))
            return newObjectHandle<NiftiImageData3DDisplacement>();
        if (boost::iequals(name, "NiftiImageData3DDeformation"))
            return newObjectHandle<NiftiImageData3DDeformation>();
        if (boost::iequals(name, "SIRFRegNiftyAladinSym"))
            return newObjectHandle<SIRFRegNiftyAladinSym<float> >();
        if (boost::iequals(name, "SIRFRegNiftyF3dSym"))
            return newObjectHandle<SIRFRegNiftyF3dSym<float> >();
        if (boost::iequals(name, "SIRFRegNiftyResample"))
            return newObjectHandle<SIRFRegNiftyResample>();
        if (boost::iequals(name, "SIRFRegImageWeightedMean"))
            return newObjectHandle<SIRFRegImageWeightedMean>();
        if (boost::iequals(name, "SIRFRegAffineTransformation"))
            return newObjectHandle<SIRFRegAffineTransformation>();
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
        if (boost::iequals(obj, "NiftiImageData"))
            return cSIRFReg_NiftiImageDataParameter(handle, name);
        if (boost::iequals(obj, "SIRFReg"))
            return cSIRFReg_SIRFRegParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegNiftyResample"))
            return cSIRFReg_SIRFRegNiftyResampleParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegImageWeightedMean"))
            return cSIRFReg_SIRFRegImageWeightedMeanParameter(handle, name);
        if (boost::iequals(obj, "SIRFRegAffineTransformation"))
            return cSIRFReg_SIRFRegAffineTransformationParameter(handle, name);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

// constructors from file
extern "C"
void* cSIRFReg_objectFromFile(const char* name, const char* filename)
{
	try {
        if (boost::iequals(name, "NiftiImageData")) {
            shared_ptr<NiftiImageData>
                sptr(new NiftiImageData(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImageData3D")) {
            shared_ptr<NiftiImageData3D>
                sptr(new NiftiImageData3D(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImageData3DTensor")) {
            shared_ptr<NiftiImageData3DTensor>
                sptr(new NiftiImageData3DTensor(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImageData3DDisplacement")) {
            shared_ptr<NiftiImageData3DDisplacement>
                sptr(new NiftiImageData3DDisplacement(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "NiftiImageData3DDeformation")) {
            shared_ptr<NiftiImageData3DDeformation>
                sptr(new NiftiImageData3DDeformation(filename));
            return newObjectHandle(sptr);
        }
        if (boost::iequals(name, "SIRFRegAffineTransformation")) {
            shared_ptr<SIRFRegAffineTransformation>
                sptr(new SIRFRegAffineTransformation(filename));
            return newObjectHandle(sptr);
        }
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImageData
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImageData_print_headers(const int num_ims, const void* im1, const void* im2, const void* im3, const void* im4, const void* im5)
{
    try {
        std::vector<NiftiImageData> vec;
        if (num_ims >= 1) vec.push_back(objectFromHandle<NiftiImageData>(im1));
        if (num_ims >= 2) vec.push_back(objectFromHandle<NiftiImageData>(im2));
        if (num_ims >= 3) vec.push_back(objectFromHandle<NiftiImageData>(im3));
        if (num_ims >= 4) vec.push_back(objectFromHandle<NiftiImageData>(im4));
        if (num_ims >= 5) vec.push_back(objectFromHandle<NiftiImageData>(im5));
        NiftiImageData::print_headers(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_save_to_file(const void* ptr, const char* filename, const int datatype)
{
	try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(ptr);
        im.save_to_file(filename,datatype);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_fill(const void* ptr, const float val)
{
    try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(ptr);
        im.fill(val);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_fill_arr(const void* ptr, size_t ptr_data)
{
    try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(ptr);
        float* data = (float*)ptr_data;
        for (int i=0; i<int(im.get_raw_nifti_sptr()->nvox); ++i)
            im(i) = data[i];
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_deep_copy(const void* copy_ptr, const void* orig_ptr)
{
    try {
        NiftiImageData& orig = objectFromHandle<NiftiImageData>(orig_ptr);
        NiftiImageData& copy = objectFromHandle<NiftiImageData>(copy_ptr);
        copy = orig.deep_copy();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_get_dimensions(const void* ptr, size_t ptr_dim)
{
    try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(ptr);
        int* dim = (int*)ptr_dim;
        for (int i=0; i<8; ++i)
            dim[i] = im.get_dimensions()[i];
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_get_data(const void* ptr, size_t ptr_data)
{
    try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(ptr);
        NiftiImageData copy = im.deep_copy();
        float* data = (float*)ptr_data;
        size_t mem = copy.get_raw_nifti_sptr()->nvox * size_t(copy.get_raw_nifti_sptr()->nbyper);
        // Copy!
        memcpy(data, copy.get_raw_nifti_sptr()->data, mem);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type)
{
    try {
        NiftiImageData& res = objectFromHandle<NiftiImageData>(res_ptr);
        NiftiImageData& im1 = objectFromHandle<NiftiImageData>(im1_ptr);
        NiftiImageData& im2 = objectFromHandle<NiftiImageData>(im2_ptr);

        if      (maths_type == 0) res = im1 + im2;
        else if (maths_type == 1) res = im1 - im2;
        else
            throw std::runtime_error("cSIRFReg_NiftiImageData_maths_im: Bad maths type (0=add, 1=subtract).");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_maths_num(const void *res_ptr, const void* im1_ptr, const float val, const int maths_type)
{
    try {
        NiftiImageData& res = objectFromHandle<NiftiImageData>(res_ptr);
        NiftiImageData& im1 = objectFromHandle<NiftiImageData>(im1_ptr);

        if      (maths_type == 0) res = im1 + val;
        else if (maths_type == 1) res = im1 - val;
        else if (maths_type == 2) res = im1 * val;
        else
            throw std::runtime_error("cSIRFReg_NiftiImageData_maths_val: Bad maths type (0=add, 1=subtract, 2=multiply.");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_equal(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImageData& im1 = objectFromHandle<NiftiImageData>(im1_ptr);
        NiftiImageData& im2 = objectFromHandle<NiftiImageData>(im2_ptr);
        return dataHandle<int>(im1 == im2);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_norm(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImageData& im1 = objectFromHandle<NiftiImageData>(im1_ptr);
        NiftiImageData& im2 = objectFromHandle<NiftiImageData>(im2_ptr);
        return dataHandle<float>(im1.get_norm(im2));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_get_original_datatype(const void* im_ptr)
{
    try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(im_ptr);
        return charDataHandleFromCharData(nifti_datatype_to_string(im.get_original_datatype()));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_crop(const void* im_ptr, size_t min_index_ptr, size_t max_index_ptr)
{
    try {
        NiftiImageData& im = objectFromHandle<NiftiImageData>(im_ptr);
        int* min_index = (int*)min_index_ptr;
        int* max_index = (int*)max_index_ptr;
        im.crop(min_index,max_index);
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      NiftiImageData3D
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImageData3D_from_PETImageData(void* ptr)
{
	try {
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(ptr);
        shared_ptr<NiftiImageData3D>
            sptr(new NiftiImageData3D(pet_im));
        return newObjectHandle(sptr);
    }
	CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3D_copy_data_to(const void* ptr, const void* obj)
{
    try {
        NiftiImageData3D& im = objectFromHandle<NiftiImageData3D>(ptr);
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(obj);
        im.copy_data_to(pet_im);
        return new DataHandle;
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      NiftiImageData3DTensor
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_save_to_file_split_xyz_components(const void *ptr, const char* filename, const int datatype)
{
	try {
        NiftiImageData3DTensor& im = objectFromHandle<NiftiImageData3DTensor>(ptr);
        im.save_to_file_split_xyz_components(filename, datatype);
		return new DataHandle;
	}
	CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_create_from_3D_image(const void* ptr, const void* obj)
{
    try {
        NiftiImageData3DTensor& im = objectFromHandle<NiftiImageData3DTensor>(ptr);
        NiftiImageData3D& im3d = objectFromHandle<NiftiImageData3D>(obj);
        im.create_from_3D_image(im3d);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr)
{
    try {
        NiftiImageData3D& x = objectFromHandle<NiftiImageData3D>(x_ptr);
        NiftiImageData3D& y = objectFromHandle<NiftiImageData3D>(y_ptr);
        NiftiImageData3D& z = objectFromHandle<NiftiImageData3D>(z_ptr);

        shared_ptr<NiftiImageData3DTensor> sptr;
        if (strcmp(obj,"NiftiImageData3DTensor"))
            sptr.reset(new NiftiImageData3DTensor(x,y,z));
        else if (strcmp(obj,"NiftiImageData3DDisplacement"))
            sptr.reset(new NiftiImageData3DDisplacement(x,y,z));
        else if (strcmp(obj,"NiftiImageData3DDeformation"))
            sptr.reset(new NiftiImageData3DDeformation(x,y,z));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_flip_component(const void *ptr, const int dim)
{
    try {
        NiftiImageData3DTensor& im = objectFromHandle<NiftiImageData3DTensor>(ptr);
        im.flip_component(dim);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImageData3DDeformation
// -------------------------------------------------------------------------------- //
void* cSIRFReg_NiftiImageData3DDeformation_compose_single_deformation(const void* im, const int num_elements, const char* types, const void* trans1, const void* trans2, const void* trans3, const void* trans4, const void* trans5)
{
    try {
        // This is an ugly hack because I can't get virtual methods to work for multiple inherited (NiftiImageData3DDeformation/NiftiImageData3DDisplacement).
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
                trans_vec.push_back(&objectFromHandle<SIRFRegAffineTransformation>(vec.at(i)));
            else if (types[i] == '2')
                trans_vec.push_back(&objectFromHandle<NiftiImageData3DDisplacement>(vec.at(i)));
            else if (types[i] == '3')
                trans_vec.push_back(&objectFromHandle<NiftiImageData3DDeformation>(vec.at(i)));

        NiftiImageData3D& ref = objectFromHandle<NiftiImageData3D>(im);
        shared_ptr<NiftiImageData3DDeformation> def_sptr
                (new NiftiImageData3DDeformation(NiftiImageData3DDeformation::compose_single_deformation(trans_vec, ref).deep_copy()));
        return newObjectHandle(def_sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DDeformation_create_from_disp(const void* ptr, const void* disp_ptr)
{
    try {
        NiftiImageData3DDeformation&  def  = objectFromHandle<NiftiImageData3DDeformation>(ptr);
        NiftiImageData3DDisplacement& disp = objectFromHandle<NiftiImageData3DDisplacement>(disp_ptr);
        def.create_from_disp(disp);
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImageData3DDisplacement
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImageData3DDisplacement_create_from_def(const void* ptr, const void* def_ptr)
{
    try {
        NiftiImageData3DDisplacement& disp = objectFromHandle<NiftiImageData3DDisplacement>(ptr);
        NiftiImageData3DDeformation&  def  = objectFromHandle<NiftiImageData3DDeformation>(def_ptr);
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
        shared_ptr<NiftiImageData3DDeformation> sptr;
        if (strcmp(transform_type, "forward_deformation") == 0)
            return newObjectHandle(shared_ptr<NiftiImageData3DDeformation>(new NiftiImageData3DDeformation(reg.get_deformation_field_forward().deep_copy())));
        else if (strcmp(transform_type, "inverse_deformation") == 0)
            return newObjectHandle(shared_ptr<NiftiImageData3DDeformation>(new NiftiImageData3DDeformation(reg.get_deformation_field_inverse().deep_copy())));
        else if (strcmp(transform_type, "forward_displacement") == 0)
            return newObjectHandle(shared_ptr<NiftiImageData3DDisplacement>(new NiftiImageData3DDisplacement(reg.get_displacement_field_forward().deep_copy())));
        else if (strcmp(transform_type, "inverse_displacement") == 0)
            return newObjectHandle(shared_ptr<NiftiImageData3DDisplacement>(new NiftiImageData3DDisplacement(reg.get_displacement_field_inverse().deep_copy())));
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
        shared_ptr<SIRFRegAffineTransformation> sptr;
        if (strcmp(dir, "forward") == 0)
            sptr.reset(new SIRFRegAffineTransformation(reg.get_transformation_matrix_forward().deep_copy()));
        else if (strcmp(dir, "inverse") == 0)
            sptr.reset(new SIRFRegAffineTransformation(reg.get_transformation_matrix_inverse().deep_copy()));
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
            res.add_transformation_affine(objectFromHandle<SIRFRegAffineTransformation>(trans));
        else if (strcmp(type, "displacement") == 0)
            res.add_transformation_disp(objectFromHandle<NiftiImageData3DDisplacement>(trans));
        else if (strcmp(type, "deformation") == 0)
            res.add_transformation_def(objectFromHandle<NiftiImageData3DDeformation>(trans));
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
        NiftiImageData& im = objectFromHandle<NiftiImageData>(obj);
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
        im_weight.add_image(NiftiImageData(filename),weight);
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

        if (strcmp(name,"SIRFRegAffineTransformation") == 0)
            trans = &objectFromHandle<SIRFRegAffineTransformation>(ptr);
        else if (strcmp(name,"NiftiImageData3DDisplacement") == 0)
            trans = &objectFromHandle<NiftiImageData3DDisplacement>(ptr);
        else if (strcmp(name,"NiftiImageData3DDeformation") == 0)
            trans = &objectFromHandle<NiftiImageData3DDeformation>(ptr);
        else
            throw std::runtime_error("cSIRFReg_SIRFRegTransformation_get_as_deformation_field: type should be affine, disp or def.");

        NiftiImageData& ref_im = objectFromHandle<NiftiImageData>(ref);
        shared_ptr<NiftiImageData3DDeformation> sptr
                (new NiftiImageData3DDeformation(trans->get_as_deformation_field(ref_im)));

        return newObjectHandle(sptr);
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SIRFRegAffineTransformation
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_construct_from_TM(size_t ptr_TM)
{
    try {
        float* TM = (float*)ptr_TM;

        SIRFRegAffineTransformation trans;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                trans[i][j] = TM[i+j*4];

        shared_ptr<SIRFRegAffineTransformation> sptr(new SIRFRegAffineTransformation(trans));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_deep_copy(const void* ptr)
{
    try {
        SIRFRegAffineTransformation& mat = objectFromHandle<SIRFRegAffineTransformation>(ptr);
        shared_ptr<SIRFRegAffineTransformation> sptr(new SIRFRegAffineTransformation(mat.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_save_to_file(const void* ptr, const char* filename)
{
    try {
        SIRFRegAffineTransformation& mat = objectFromHandle<SIRFRegAffineTransformation>(ptr);
        mat.save_to_file(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_as_array(const void* ptr, size_t ptr_TM)
{
    try {
        SIRFRegAffineTransformation& tm = objectFromHandle<SIRFRegAffineTransformation>(ptr);
        float* TM = (float*)ptr_TM;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                TM[i*4+j] = tm[i][j];
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_get_identity()
{
    try {
        shared_ptr<SIRFRegAffineTransformation> sptr(new SIRFRegAffineTransformation(SIRFRegAffineTransformation::get_identity()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_get_inverse(const void* ptr)
{
    try {
        SIRFRegAffineTransformation& tm = objectFromHandle<SIRFRegAffineTransformation>(ptr);
        shared_ptr<SIRFRegAffineTransformation> sptr(new SIRFRegAffineTransformation(tm.get_inverse()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_mul(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        SIRFRegAffineTransformation& mat1 = objectFromHandle<SIRFRegAffineTransformation>(mat1_ptr);
        SIRFRegAffineTransformation& mat2 = objectFromHandle<SIRFRegAffineTransformation>(mat2_ptr);
        shared_ptr<SIRFRegAffineTransformation> sptr(new SIRFRegAffineTransformation(mat1*mat2));
        return newObjectHandle(sptr);
    }
    CATCH;
}

extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_equal(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        SIRFRegAffineTransformation& mat1 = objectFromHandle<SIRFRegAffineTransformation>(mat1_ptr);
        SIRFRegAffineTransformation& mat2 = objectFromHandle<SIRFRegAffineTransformation>(mat2_ptr);
        return dataHandle<int>(mat1 == mat2);
    }
    CATCH;
}
