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

#include "sirf/iUtilities/DataHandle.h"
#include "cReg.h"
#include "sirf/cReg/csirfreg_p.h"
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/NiftiImageData3DTensor.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"
#include "sirf/cReg/NiftyAladinSym.h"
#include "sirf/cReg/NiftyF3dSym.h"
#include "sirf/cReg/NiftyResample.h"
#include "sirf/cReg/ImageWeightedMean.h"
#include "sirf/cReg/Transformation.h"
#include "sirf/cReg/AffineTransformation.h"

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
        if (strcmp(name, "NiftiImageData") == 0)
            return newObjectHandle(std::shared_ptr<NiftiImageData<float> >(new NiftiImageData<float>));
        if (strcmp(name, "NiftiImageData3D") == 0)
            return newObjectHandle(std::shared_ptr<NiftiImageData3D<float> >(new NiftiImageData3D<float>));
        if (strcmp(name, "NiftiImageData3DTensor") == 0)
            return newObjectHandle(std::shared_ptr<NiftiImageData3DTensor<float> >(new NiftiImageData3DTensor<float>));
        if (strcmp(name, "NiftiImageData3DDisplacement") == 0)
            return newObjectHandle(std::shared_ptr<NiftiImageData3DDisplacement<float> >(new NiftiImageData3DDisplacement<float>));
        if (strcmp(name, "NiftiImageData3DDeformation") == 0)
            return newObjectHandle(std::shared_ptr<NiftiImageData3DDeformation<float> >(new NiftiImageData3DDeformation<float>));
        if (strcmp(name, "SIRFRegNiftyAladinSym") == 0)
            return newObjectHandle(std::shared_ptr<NiftyAladinSym<float> >(new NiftyAladinSym<float>));
        if (strcmp(name, "SIRFRegNiftyF3dSym") == 0)
            return newObjectHandle(std::shared_ptr<NiftyF3dSym<float> >(new NiftyF3dSym<float>));
        if (strcmp(name, "SIRFRegNiftyResample") == 0)
            return newObjectHandle(std::shared_ptr<NiftyResample<float> >(new NiftyResample<float>));
        if (strcmp(name, "SIRFRegImageWeightedMean") == 0)
            return newObjectHandle(std::shared_ptr<ImageWeightedMean<float> >(new ImageWeightedMean<float>));
        if (strcmp(name, "SIRFRegAffineTransformation") == 0)
            return newObjectHandle(std::shared_ptr<AffineTransformation<float> >(new AffineTransformation<float>));
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
        if (strcmp(obj, "SIRFReg") == 0)
            return cSIRFReg_setSIRFRegParameter(ptr_s, name, ptr_v);
        if (strcmp(obj, "SIRFRegNiftyF3dSym") == 0)
            return cSIRFReg_setSIRFRegNiftyF3dSymParameter(ptr_s, name, ptr_v);
        if (strcmp(obj, "SIRFRegNiftyResample") == 0)
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
        if (strcmp(obj, "NiftiImageData") == 0)
            return cSIRFReg_NiftiImageDataParameter(handle, name);
        if (strcmp(obj, "SIRFReg") == 0)
            return cSIRFReg_SIRFRegParameter(handle, name);
        if (strcmp(obj, "SIRFRegNiftyResample") == 0)
            return cSIRFReg_SIRFRegNiftyResampleParameter(handle, name);
        if (strcmp(obj, "SIRFRegImageWeightedMean") == 0)
            return cSIRFReg_SIRFRegImageWeightedMeanParameter(handle, name);
        if (strcmp(obj, "SIRFRegAffineTransformation") == 0)
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
        if (strcmp(name, "NiftiImageData") == 0) {
            std::shared_ptr<NiftiImageData<float> >
                sptr(new NiftiImageData<float>(filename));
            return newObjectHandle(sptr);
        }
        if (strcmp(name, "NiftiImageData3D") == 0) {
            std::shared_ptr<NiftiImageData3D<float> >
                sptr(new NiftiImageData3D<float>(filename));
            return newObjectHandle(sptr);
        }
        if (strcmp(name, "NiftiImageData3DTensor") == 0) {
            std::shared_ptr<NiftiImageData3DTensor<float> >
                sptr(new NiftiImageData3DTensor<float>(filename));
            return newObjectHandle(sptr);
        }
        if (strcmp(name, "NiftiImageData3DDisplacement") == 0) {
            std::shared_ptr<NiftiImageData3DDisplacement<float> >
                sptr(new NiftiImageData3DDisplacement<float>(filename));
            return newObjectHandle(sptr);
        }
        if (strcmp(name, "NiftiImageData3DDeformation") == 0) {
            std::shared_ptr<NiftiImageData3DDeformation<float> >
                sptr(new NiftiImageData3DDeformation<float>(filename));
            return newObjectHandle(sptr);
        }
        if (strcmp(name, "SIRFRegAffineTransformation") == 0) {
            std::shared_ptr<AffineTransformation<float> >
                sptr(new AffineTransformation<float>(filename));
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
        std::vector<NiftiImageData<float> > vec;
        if (num_ims >= 1) vec.push_back(objectFromHandle<NiftiImageData<float> >(im1));
        if (num_ims >= 2) vec.push_back(objectFromHandle<NiftiImageData<float> >(im2));
        if (num_ims >= 3) vec.push_back(objectFromHandle<NiftiImageData<float> >(im3));
        if (num_ims >= 4) vec.push_back(objectFromHandle<NiftiImageData<float> >(im4));
        if (num_ims >= 5) vec.push_back(objectFromHandle<NiftiImageData<float> >(im5));
        std::vector<const NiftiImageData<float>*> vec_ptr;
        for (int i=0; i<vec.size(); ++i)
            vec_ptr.push_back(&vec[i]);
        NiftiImageData<float>::print_headers(vec_ptr);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_write(const void* ptr, const char* filename, const int datatype)
{
	try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(ptr);
        im.write(filename,datatype);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_fill(const void* ptr, const float val)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(ptr);
        im.fill(val);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_fill_arr(const void* ptr, size_t ptr_data)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(ptr);
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
        NiftiImageData<float>& orig = objectFromHandle<NiftiImageData<float> >(orig_ptr);
        NiftiImageData<float>& copy = objectFromHandle<NiftiImageData<float> >(copy_ptr);
        copy = orig;
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_get_dimensions(const void* ptr, size_t ptr_dim)
{
    try {
        NiftiImageData<float> & im = objectFromHandle<NiftiImageData<float> >(ptr);
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
        const NiftiImageData<float>& im = objectFromHandle<const NiftiImageData<float> >(ptr);
        float* data = (float*)ptr_data;
        size_t mem = im.get_raw_nifti_sptr()->nvox * size_t(im.get_raw_nifti_sptr()->nbyper);
        // Copy!
        memcpy(data, im.get_raw_nifti_sptr()->data, mem);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type)
{
    try {
        NiftiImageData<float>& res = objectFromHandle<NiftiImageData<float> >(res_ptr);
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);

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
        NiftiImageData<float>& res = objectFromHandle<NiftiImageData<float> >(res_ptr);
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);

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
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);
        return dataHandle<int>(im1 == im2);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_norm(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);
        return dataHandle<float>(im1.get_norm(im2));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_get_original_datatype(const void* im_ptr)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(im_ptr);
        return charDataHandleFromCharData(nifti_datatype_to_string(im.get_original_datatype()));
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData_crop(const void* im_ptr, size_t min_index_ptr, size_t max_index_ptr)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(im_ptr);
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
/* TODO UNCOMMENT WHEN GEOMETRICAL INFO IS IMPLEMENTED
extern "C"
void* cSIRFReg_NiftiImageData3D_from_PETImageData(void* ptr)
{
	try {
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(ptr);
        shared_ptr<NiftiImageData3D<float> >
            sptr(new NiftiImageData3D<float>(pet_im));
        return newObjectHandle(sptr);
    }
	CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3D_copy_data_to(const void* ptr, const void* obj)
{
    try {
        NiftiImageData3D<float>& im = objectFromHandle<NiftiImageData3D<float> >(ptr);
        sirf::PETImageData& pet_im = objectFromHandle<sirf::PETImageData>(obj);
        im.copy_data_to(pet_im);
        return new DataHandle;
    }
    CATCH;
}
*/
// -------------------------------------------------------------------------------- //
//      NiftiImageData3DTensor
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_write_split_xyz_components(const void *ptr, const char* filename, const int datatype)
{
	try {
        NiftiImageData3DTensor<float>& im = objectFromHandle<NiftiImageData3DTensor<float> >(ptr);
        im.write_split_xyz_components(filename, datatype);
		return new DataHandle;
	}
	CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_create_from_3D_image(const void* ptr, const void* obj)
{
    try {
        NiftiImageData3DTensor<float>& im = objectFromHandle<NiftiImageData3DTensor<float> >(ptr);
        NiftiImageData3D<float>& im3d = objectFromHandle<NiftiImageData3D<float> >(obj);
        im.create_from_3D_image(im3d);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr)
{
    try {
        NiftiImageData3D<float>& x = objectFromHandle<NiftiImageData3D<float> >(x_ptr);
        NiftiImageData3D<float>& y = objectFromHandle<NiftiImageData3D<float> >(y_ptr);
        NiftiImageData3D<float>& z = objectFromHandle<NiftiImageData3D<float> >(z_ptr);

        std::shared_ptr<NiftiImageData3DTensor<float> > sptr;
        if (strcmp(obj,"NiftiImageData3DTensor") == 0)
            sptr.reset(new NiftiImageData3DTensor<float>(x,y,z));
        else if (strcmp(obj,"NiftiImageData3DDisplacement") == 0)
            sptr.reset(new NiftiImageData3DDisplacement<float>(x,y,z));
        else if (strcmp(obj,"NiftiImageData3DDeformation") == 0)
            sptr.reset(new NiftiImageData3DDeformation<float>(x,y,z));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DTensor_flip_component(const void *ptr, const int dim)
{
    try {
        NiftiImageData3DTensor<float>& im = objectFromHandle<NiftiImageData3DTensor<float> >(ptr);
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
        std::vector<const Transformation<float> *> trans_vec;
        for (int i=0; i<num_elements; ++i)
            if      (types[i] == '1')
                trans_vec.push_back(&objectFromHandle<const AffineTransformation<float> >(vec.at(i)));
            else if (types[i] == '2')
                trans_vec.push_back(&objectFromHandle<const NiftiImageData3DDisplacement<float> >(vec.at(i)));
            else if (types[i] == '3')
                trans_vec.push_back(&objectFromHandle<const NiftiImageData3DDeformation<float> >(vec.at(i)));

        const NiftiImageData3D<float>& ref = objectFromHandle<const NiftiImageData3D<float> >(im);
        const std::shared_ptr<const NiftiImageData3DDeformation<float> > def_sptr
                (new const NiftiImageData3DDeformation<float>(NiftiImageData3DDeformation<float>::compose_single_deformation(trans_vec, ref)));
        return newObjectHandle(def_sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_NiftiImageData3DDeformation_create_from_disp(const void* disp_ptr)
{
    try {
        NiftiImageData3DDisplacement<float>& disp = objectFromHandle<NiftiImageData3DDisplacement<float> >(disp_ptr);
        return newObjectHandle(std::make_shared<NiftiImageData3DDeformation<float> >(disp));
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImageData3DDisplacement
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_NiftiImageData3DDisplacement_create_from_def(const void* def_ptr)
{
    try {
        NiftiImageData3DDeformation<float>& def = objectFromHandle<NiftiImageData3DDeformation<float> >(def_ptr);
        return newObjectHandle(std::make_shared<NiftiImageData3DDisplacement<float> >(def));
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
        Registration<float>& reg = objectFromHandle<Registration<float> >(ptr);
        reg.process();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_get_deformation_displacement_image(const void* ptr, const char *transform_type)
{
    try {
        Registration<float>& reg = objectFromHandle<Registration<float>>(ptr);
        if (strcmp(transform_type, "forward_deformation") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> >(reg.get_deformation_field_forward()));
        else if (strcmp(transform_type, "inverse_deformation") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> >(reg.get_deformation_field_inverse()));
        else if (strcmp(transform_type, "forward_displacement") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(reg.get_displacement_field_forward()));
        else if (strcmp(transform_type, "inverse_displacement") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(reg.get_displacement_field_inverse()));
        else
            throw std::runtime_error("cSIRFReg_SIRFReg_get_deformation_displacement_image: Bad return type.");
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFReg_set_parameter(const void* ptr, const char* par, const char* arg1, const char* arg2)
{
    try {
        Registration<float>& reg = objectFromHandle<Registration<float> >(ptr);
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
        NiftyAladinSym<float>& reg = objectFromHandle<NiftyAladinSym<float> >(ptr);
        std::shared_ptr<const AffineTransformation<float> > sptr;
        if (strcmp(dir, "forward") == 0)
            sptr = reg.get_transformation_matrix_forward();
        else if (strcmp(dir, "inverse") == 0)
            sptr = reg.get_transformation_matrix_inverse();
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
        NiftyResample<float>& res = objectFromHandle<NiftyResample<float> >(self);
        if (strcmp(type, "affine") == 0)
            res.add_transformation(std::make_shared<const AffineTransformation<float> >(objectFromHandle<AffineTransformation<float> >(trans)));
        else if (strcmp(type, "displacement") == 0)
            res.add_transformation(std::make_shared<const NiftiImageData3DDisplacement<float> >(objectFromHandle<NiftiImageData3DDisplacement<float> >(trans)));
        else if (strcmp(type, "deformation") == 0)
            res.add_transformation(std::make_shared<const NiftiImageData3DDeformation<float> >(objectFromHandle<NiftiImageData3DDeformation<float> >(trans)));
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
        NiftyResample<float>& res = objectFromHandle<NiftyResample<float> >(ptr);
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
        ImageWeightedMean<float>& im_weight = objectFromHandle<ImageWeightedMean<float> >(ptr);
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(obj);
        im_weight.add_image(im,weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight)
{
    try {
        ImageWeightedMean<float>& im_weight = objectFromHandle<ImageWeightedMean<float> >(ptr);
        im_weight.add_image(NiftiImageData<float>(filename),weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegImageWeightedMean_process(void* ptr)
{
    try {
        ImageWeightedMean<float>& im_weight = objectFromHandle<ImageWeightedMean<float> >(ptr);
        im_weight.process();
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      Transformation
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegTransformation_get_as_deformation_field(const void* ptr, const char* name, const void* ref)
{
    try {
        Transformation<float> *trans;

        if (strcmp(name,"SIRFRegAffineTransformation") == 0)
            trans = &objectFromHandle<AffineTransformation<float> >(ptr);
        else if (strcmp(name,"NiftiImageData3DDisplacement") == 0)
            trans = &objectFromHandle<NiftiImageData3DDisplacement<float> >(ptr);
        else if (strcmp(name,"NiftiImageData3DDeformation") == 0)
            trans = &objectFromHandle<NiftiImageData3DDeformation<float> >(ptr);
        else
            throw std::runtime_error("cSIRFReg_Transformation_get_as_deformation_field: type should be affine, disp or def.");

        NiftiImageData<float>& ref_im = objectFromHandle<NiftiImageData<float> >(ref);
        std::shared_ptr<NiftiImageData3DDeformation<float> > sptr
                (new NiftiImageData3DDeformation<float>(trans->get_as_deformation_field(ref_im)));

        return newObjectHandle(sptr);
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      AffineTransformation
// -------------------------------------------------------------------------------- //
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_construct_from_TM(size_t ptr_TM)
{
    try {
        float* TM = (float*)ptr_TM;

        AffineTransformation<float> trans;
        for (int i=0; i<4; ++i)
            for (int j=0; j<4; ++j)
                trans[i][j] = TM[i+j*4];

        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>(trans));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_deep_copy(const void* ptr)
{
    try {
        AffineTransformation<float>& mat = objectFromHandle<AffineTransformation<float> >(ptr);
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>(mat.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_write(const void* ptr, const char* filename)
{
    try {
        AffineTransformation<float>& mat = objectFromHandle<AffineTransformation<float> >(ptr);
        mat.write(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_as_array(const void* ptr, size_t ptr_TM)
{
    try {
        AffineTransformation<float>& tm = objectFromHandle<AffineTransformation<float> >(ptr);
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
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>);
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_get_inverse(const void* ptr)
{
    try {
        AffineTransformation<float>& tm = objectFromHandle<AffineTransformation<float> >(ptr);
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>(tm.get_inverse()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_mul(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        AffineTransformation<float>& mat1 = objectFromHandle<AffineTransformation<float> >(mat1_ptr);
        AffineTransformation<float>& mat2 = objectFromHandle<AffineTransformation<float> >(mat2_ptr);
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>(mat1*mat2));
        return newObjectHandle(sptr);
    }
    CATCH;
}

extern "C"
void* cSIRFReg_SIRFRegAffineTransformation_equal(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        AffineTransformation<float>& mat1 = objectFromHandle<AffineTransformation<float> >(mat1_ptr);
        AffineTransformation<float>& mat2 = objectFromHandle<AffineTransformation<float> >(mat2_ptr);
        return dataHandle<int>(mat1 == mat2);
    }
    CATCH;
}
