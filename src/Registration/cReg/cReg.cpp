/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2017 - 2020 University College London

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
#include "sirf/Reg/cReg.h"
#include "sirf/Reg/cReg_p.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftyAladinSym.h"
#include "sirf/Reg/NiftyF3dSym.h"
#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/ImageWeightedMean.h"
#include "sirf/Reg/Transformation.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/Quaternion.h"
#ifdef SIRF_SPM
#include "sirf/Reg/SPMRegistration.h"
#endif

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
void* cReg_newObject(const char* name)
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
        if (strcmp(name, "NiftyAladinSym") == 0)
            return newObjectHandle(std::shared_ptr<NiftyAladinSym<float> >(new NiftyAladinSym<float>));
        if (strcmp(name, "NiftyF3dSym") == 0)
            return newObjectHandle(std::shared_ptr<NiftyF3dSym<float> >(new NiftyF3dSym<float>));
        if (strcmp(name, "NiftyResample") == 0)
            return newObjectHandle(std::shared_ptr<NiftyResample<float> >(new NiftyResample<float>));
        if (strcmp(name, "ImageWeightedMean") == 0)
            return newObjectHandle(std::shared_ptr<ImageWeightedMean<float> >(new ImageWeightedMean<float>));
        if (strcmp(name, "AffineTransformation") == 0)
            return newObjectHandle(std::shared_ptr<AffineTransformation<float> >(new AffineTransformation<float>));
#ifdef SIRF_SPM
        if (strcmp(name, "SPMRegistration") == 0)
            return newObjectHandle(std::shared_ptr<SPMRegistration<float> >(new SPMRegistration<float>));
#endif
		return unknownObject("object", name, __FILE__, __LINE__);
	}
	CATCH;
}

// set parameters
extern "C"
void* setParameter
(void* ptr_s, const char* obj, const char* name, const void* ptr_v)
{
	try {
        if (strcmp(obj, "Registration") == 0)
            return cReg_setRegistrationParameter(ptr_s, name, ptr_v);
        if (strcmp(obj, "NiftyRegistration") == 0)
            return cReg_setNiftyRegistrationParameter(ptr_s, name, ptr_v);
#ifdef SIRF_SPM
        if (strcmp(obj, "SPMRegistration") == 0)
            return cReg_setSPMRegistrationParameter(ptr_s, name, ptr_v);
#endif
        if (strcmp(obj, "NiftyF3dSym") == 0)
            return cReg_setNiftyF3dSymParameter(ptr_s, name, ptr_v);
        if (strcmp(obj, "NiftyResample") == 0)
            return cReg_setNiftyResampleParameter(ptr_s, name, ptr_v);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

// get functions
extern "C"
void* parameter(const void* ptr, const char* obj, const char* name)
{
	try {
		CAST_PTR(DataHandle, handle, ptr);
        if (strcmp(obj, "NiftiImageData") == 0)
            return cReg_NiftiImageDataParameter(handle, name);
        if (strcmp(obj, "NiftyResample") == 0)
            return cReg_NiftyResampleParameter(handle, name);
        if (strcmp(obj, "ImageWeightedMean") == 0)
            return cReg_ImageWeightedMeanParameter(handle, name);
        if (strcmp(obj, "AffineTransformation") == 0)
            return cReg_AffineTransformationParameter(handle, name);
		return unknownObject("object", obj, __FILE__, __LINE__);
	}
	CATCH;
}

// constructors from file
extern "C"
void* cReg_objectFromFile(const char* name, const char* filename)
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
        if (strcmp(name, "AffineTransformation") == 0) {
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
void* cReg_NiftiImageData_print_headers(const void* handle_vector_ptr)
{
    try {
        const DataHandleVector handle_vector = objectFromHandle<const DataHandleVector>(handle_vector_ptr);
        std::vector<const NiftiImageData<float>*> vec;
        for (unsigned i=0; i<handle_vector.size(); ++i)
            vec.push_back(&objectFromHandle<const NiftiImageData<float> >(handle_vector.at(i)));
        NiftiImageData<float>::print_headers(vec);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_write(const void* ptr, const char* filename, const int datatype)
{
	try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(ptr);
        im.write(filename,datatype);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_fill(const void* ptr, const float val)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(ptr);
        im.fill(val);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_fill_arr(const void* ptr, size_t ptr_data)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(ptr);

        const int *dims = im.get_dimensions();
        int dim_x = dims[1];
        int dim_y = dims[2];
        int dim_z = dims[3];
        int dim_t = dims[4];
        int dim_u = dims[5];
        int dim_v = dims[6];
        int dim_w = dims[7];
        // Only implemented for 3D scalar or tensor images. of x,y,z,t,u,v,w, throw an error if t,v,w are != 1.
        if (dim_t!=1 || dim_v!=1 || dim_w!=1)
            throw std::runtime_error("fill only implemented for 3D scalar or tensor images (should be easy to extend).");

        // Get arrays
        float *im_data = static_cast<float*>(im.get_raw_nifti_sptr()->data);
        float *data = (float*)ptr_data;

        // nifti_image data are stored as u,x,y,z, whereas python and matlab need x,y,z,u
        int wrap_idx;
        for (int u=0; u<dim_u; ++u) {
            for (int x=0; x<dim_x; ++x) {
                for (int y=0; y<dim_y; ++y) {
                    for (int z=0; z<dim_z; ++z) {
                        int nifti_idx[7] = { x,y,z,0,u,0,0 };
                        wrap_idx  = u + dim_u*(x + dim_x*(y + dim_y*(z)));
                        im(nifti_idx) = data[wrap_idx];
                    }
                }
            }
        }

        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_deep_copy(const void* copy_ptr, const void* orig_ptr)
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
void* cReg_NiftiImageData_get_dimensions(const void* ptr, size_t ptr_dim)
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
void* cReg_NiftiImageData_get_voxel_sizes(const void* ptr, PTR_FLOAT ptr_out)
{
    try {
        NiftiImageData<float> & im = objectFromHandle<NiftiImageData<float> >(ptr);
        float* dim = (float*)ptr_out;
        for (int i=0; i<8; ++i)
            dim[i] = im.get_raw_nifti_sptr()->pixdim[i];
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_as_array(const void* ptr, size_t ptr_data)
{
    try {
        const NiftiImageData<float>& im = objectFromHandle<const NiftiImageData<float> >(ptr);
        const int *dims = im.get_dimensions();
        int dim_x = dims[1];
        int dim_y = dims[2];
        int dim_z = dims[3];
        int dim_t = dims[4];
        int dim_u = dims[5];
        int dim_v = dims[6];
        int dim_w = dims[7];
        // Only implemented for 3D scalar or tensor images. of x,y,z,t,u,v,w, throw an error if t,v,w are != 1.
        if (dim_t!=1 || dim_v!=1 || dim_w!=1)
            throw std::runtime_error("as_array only implemented for 3D scalar or tensor images (should be easy to extend).");

        // Get arrays
        const float *im_data = static_cast<float*>(im.get_raw_nifti_sptr()->data);
        float *data = (float*)ptr_data;

        // nifti_image data are stored as u,x,y,z, whereas python and matlab need x,y,z,u
        int wrap_idx;
        for (int u=0; u<dim_u; ++u) {
            for (int x=0; x<dim_x; ++x) {
                for (int y=0; y<dim_y; ++y) {
                    for (int z=0; z<dim_z; ++z) {
                        int nifti_idx[7] = { x,y,z,0,u,0,0 };
                        wrap_idx  = x + dim_x*(y + dim_y*(z + dim_z*(u)));
                        data[wrap_idx] = im(nifti_idx);
                    }
                }
            }
        }
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_maths_im(const void *res_ptr, const void* im1_ptr, const void* im2_ptr, const int maths_type)
{
    try {
        NiftiImageData<float>& res = objectFromHandle<NiftiImageData<float> >(res_ptr);
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);

        if      (maths_type == 0) res = im1 + im2;
        else if (maths_type == 1) res = im1 - im2;
        else
            throw std::runtime_error("cReg_NiftiImageData_maths_im: Bad maths type (0=add, 1=subtract).");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_maths_num(const void *res_ptr, const void* im1_ptr, const float val, const int maths_type)
{
    try {
        NiftiImageData<float>& res = objectFromHandle<NiftiImageData<float> >(res_ptr);
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);

        if      (maths_type == 0) res = im1 + val;
        else if (maths_type == 1) res = im1 - val;
        else if (maths_type == 2) res = im1 * val;
        else
            throw std::runtime_error("cReg_NiftiImageData_maths_val: Bad maths type (0=add, 1=subtract, 2=multiply.");
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_equal(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);
        return dataHandle<int>(im1 == im2);
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_norm(const void* im1_ptr, const void* im2_ptr)
{
    try {
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);
        return dataHandle<float>(im1.get_norm(im2));
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_get_original_datatype(const void* im_ptr)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(im_ptr);
        return dataHandle<int>(im.get_original_datatype());
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_crop(const void* im_ptr, size_t min_index_ptr, size_t max_index_ptr)
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

extern "C"
void* cReg_NiftiImageData_set_voxel_spacing(const void* im_ptr, const float x, const float y, const float z, const int interpolation_order)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(im_ptr);
        float spacing[3] = {x,y,z};
        im.set_voxel_spacing(spacing,interpolation_order);
        return new DataHandle;
    }
    CATCH;
}

extern "C"
void* cReg_NiftiImageData_normalise_zero_and_one(const void* im_ptr)
{
    try {
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(im_ptr);
        im.normalise_zero_and_one();
        return new DataHandle;
    }
    CATCH;
}

extern "C"
void* cReg_NiftiImageData_standardise(const void* im_ptr)
{
    try{
        NiftiImageData<float>& im = objectFromHandle<NiftiImageData<float> >(im_ptr);
        im.standardise();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_get_inner_product(const void* im1_ptr, const void* im2_ptr)
{
    try{
        NiftiImageData<float>& im1 = objectFromHandle<NiftiImageData<float> >(im1_ptr);
        NiftiImageData<float>& im2 = objectFromHandle<NiftiImageData<float> >(im2_ptr);
        return dataHandle<float>(im1.get_inner_product(im2));
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData_from_SIRFImageData(void* ptr)
{
	try {
        ImageData& sirf_im = objectFromHandle<ImageData>(ptr);
        std::shared_ptr<NiftiImageData<float> >
            sptr(new NiftiImageData<float>(sirf_im));
        return newObjectHandle(sptr);
    }
	CATCH;
}

extern "C"
void* cReg_NiftiImageData_from_complex_ImageData_real_component(void* in_ptr)
{
    try {
        std::shared_ptr<ImageData> in_sptr;
        getObjectSptrFromHandle<ImageData>(in_ptr, in_sptr);
        std::shared_ptr<NiftiImageData<float> > out_sptr;
        NiftiImageData<float>::construct_NiftiImageData_from_complex_im_real_component(out_sptr, in_sptr);
        return newObjectHandle(out_sptr);
    }
	CATCH;
}

extern "C"
void* cReg_NiftiImageData_from_complex_ImageData_imag_component(void* in_ptr)
{
    try {
        std::shared_ptr<ImageData> in_sptr;
        getObjectSptrFromHandle<ImageData>(in_ptr, in_sptr);
        std::shared_ptr<NiftiImageData<float> > out_sptr;
        NiftiImageData<float>::construct_NiftiImageData_from_complex_im_imag_component(out_sptr, in_sptr);
        return newObjectHandle(out_sptr);
    }
	CATCH;
}

extern "C"
void* cReg_NiftiImageData_are_equal_to_given_accuracy(void* im1_ptr, void* im2_ptr, const float accuracy)
{
    try {
        std::shared_ptr<NiftiImageData<float> > im1_sptr, im2_sptr;
        getObjectSptrFromHandle<NiftiImageData<float> >(im1_ptr, im1_sptr);
        getObjectSptrFromHandle<NiftiImageData<float> >(im2_ptr, im2_sptr);
        return dataHandle<int>(NiftiImageData<float>::are_equal_to_given_accuracy(im1_sptr, im2_sptr, accuracy));
    }
	CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImageData3DTensor
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_NiftiImageData3DTensor_write_split_xyz_components(const void *ptr, const char* filename, const int datatype)
{
	try {
        NiftiImageData3DTensor<float>& im = objectFromHandle<NiftiImageData3DTensor<float> >(ptr);
        im.write_split_xyz_components(filename, datatype);
		return new DataHandle;
	}
	CATCH;
}
extern "C"
void* cReg_NiftiImageData3DTensor_create_from_3D_image(const void* ptr, const void* obj)
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
void* cReg_NiftiImageData3DTensor_construct_from_3_components(const char* obj, const void *x_ptr, const void *y_ptr, const void *z_ptr)
{
    try {
        ImageData& x = objectFromHandle<ImageData>(x_ptr);
        ImageData& y = objectFromHandle<ImageData>(y_ptr);
        ImageData& z = objectFromHandle<ImageData>(z_ptr);

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
void* cReg_NiftiImageData3DTensor_flip_component(const void *ptr, const int dim)
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
void* cReg_NiftiImageData3DDeformation_compose_single_deformation(const void* im, const char* types, const void* trans_vector_ptr)
{
    try {
        // This is an ugly hack because I can't get virtual methods to work for multiple inherited (NiftiImageData3DDeformation/NiftiImageData3DDisplacement).
        // So, we also give a string which tells us what type they are, and we change the template type of objectFromHandle accordingly.

        // Sorry this is so ugly.

        // There's always going to be at least two transformations, so start by putting them in the vector
        const DataHandleVector vec = objectFromHandle<const DataHandleVector>(trans_vector_ptr);
        std::vector<const Transformation<float> *> trans_vec;
        for (unsigned i=0; i<vec.size(); ++i)
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
void* cReg_NiftiImageData3DDeformation_create_from_disp(const void* disp_ptr)
{
    try {
        NiftiImageData3DDisplacement<float>& disp = objectFromHandle<NiftiImageData3DDisplacement<float> >(disp_ptr);
        return newObjectHandle(std::make_shared<NiftiImageData3DDeformation<float> >(disp));
    }
    CATCH;
}
extern "C"
void* cReg_NiftiImageData3DDeformation_get_inverse(const void* def_ptr, const void* floating_ptr)
{
    try {
        NiftiImageData3DDeformation<float>& def = objectFromHandle<NiftiImageData3DDeformation<float> >(def_ptr);
        std::shared_ptr<const NiftiImageData<float> > flo_sptr;
        getObjectSptrFromHandle<const NiftiImageData<float> >(floating_ptr, flo_sptr);
        return newObjectHandle(def.get_inverse(flo_sptr));
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftiImageData3DDisplacement
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_NiftiImageData3DDisplacement_create_from_def(const void* def_ptr)
{
    try {
        NiftiImageData3DDeformation<float>& def = objectFromHandle<NiftiImageData3DDeformation<float> >(def_ptr);
        return newObjectHandle(std::make_shared<NiftiImageData3DDisplacement<float> >(def));
    }
    CATCH;
}

// -------------------------------------------------------------------------------- //
//      Registration
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_Registration_process(void* ptr)
{
    try {
        Registration<float>& reg = objectFromHandle<Registration<float> >(ptr);
        reg.process();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_Registration_get_deformation_displacement_image(const void* ptr, const char *transform_type, const int idx)
{
    try {
        Registration<float>& reg = objectFromHandle<Registration<float>>(ptr);
        if (strcmp(transform_type, "forward_deformation") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> >(reg.get_deformation_field_forward_sptr(unsigned(idx))));
        else if (strcmp(transform_type, "inverse_deformation") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDeformation<float> >(reg.get_deformation_field_inverse_sptr(unsigned(idx))));
        else if (strcmp(transform_type, "forward_displacement") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(reg.get_displacement_field_forward_sptr(unsigned(idx))));
        else if (strcmp(transform_type, "inverse_displacement") == 0)
            return newObjectHandle(std::dynamic_pointer_cast<const NiftiImageData3DDisplacement<float> >(reg.get_displacement_field_inverse_sptr(unsigned(idx))));
        else
            throw std::runtime_error("cReg_Registration_get_deformation_displacement_image: Bad return type.");
    }
    CATCH;
}
extern "C"
void* cReg_Registration_add_floating(const void* ptr, const void* im_ptr)
{
    Registration<float>& reg = objectFromHandle<Registration<float> >(ptr);
    std::shared_ptr<const ImageData> im_sptr;
    getObjectSptrFromHandle<const ImageData>(im_ptr, im_sptr);
    reg.add_floating_image(im_sptr);
    return new DataHandle;
}
extern "C"
void* cReg_Registration_clear_floatings(const void* ptr)
{
    Registration<float>& reg = objectFromHandle<Registration<float> >(ptr);
    reg.clear_floating_images();
    return new DataHandle;
}
extern "C"
void* cReg_Registration_get_output(const void* ptr,const int idx)
{
    Registration<float>& reg = objectFromHandle<Registration<float> >(ptr);
    return newObjectHandle(reg.get_output_sptr(unsigned(idx)));
}
// -------------------------------------------------------------------------------- //
//      NiftyRegistration
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_NiftyRegistration_set_parameter(const void* ptr, const char* par, const char* arg1, const char* arg2)
{
    try {
        NiftyRegistration<float>& reg = objectFromHandle<NiftyRegistration<float> >(ptr);
        reg.set_parameter(par, arg1, arg2);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftyRegistration_print_all_wrapped_methods(const char* name)
{
    try {
        if (strcmp(name, "NiftyAladinSym") == 0)
            NiftyAladinSym<float>::print_all_wrapped_methods();
        else if (strcmp(name, "NiftyF3dSym") == 0)
            NiftyF3dSym<float>::print_all_wrapped_methods();
        else
            throw std::runtime_error("cReg_Registration_print_all_wrapped_methods: Non-existent reconstruction algorithm name.");
        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      NiftyAladinSym
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_NiftyAladin_get_TM(const void* ptr, const char* dir)
{
    try {
        NiftyAladinSym<float>& reg = objectFromHandle<NiftyAladinSym<float> >(ptr);
        std::shared_ptr<const AffineTransformation<float> > sptr;
        if (strcmp(dir, "forward") == 0)
            sptr = reg.get_transformation_matrix_forward_sptr();
        else if (strcmp(dir, "inverse") == 0)
            sptr = reg.get_transformation_matrix_inverse_sptr();
        else
            throw std::runtime_error("only accept forward or inverse as argument to dir for saving transformation matrix");
        return newObjectHandle(sptr);
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      SPM
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_SPMRegistration_get_TM(const void* ptr, const char* dir, const int idx)
{
#ifdef SIRF_SPM
    try {
        SPMRegistration<float>& reg = objectFromHandle<SPMRegistration<float> >(ptr);
        std::shared_ptr<const AffineTransformation<float> > sptr;
        if (strcmp(dir, "forward") == 0)
            sptr = reg.get_transformation_matrix_forward_sptr(unsigned(idx));
        else if (strcmp(dir, "inverse") == 0)
            sptr = reg.get_transformation_matrix_inverse_sptr(unsigned(idx));
        else
            throw std::runtime_error("only accept forward or inverse as argument to dir for saving transformation matrix");
        return newObjectHandle(sptr);
    }
    CATCH;
#else
    throw std::runtime_error("cReg_SPMRegistration_get_TM: SPM not present, you shouldn't be here.");
#endif
}
// -------------------------------------------------------------------------------- //
//      NiftyResample
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_NiftyResample_add_transformation(void* self, const void* trans, const char *type)
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
void* cReg_NiftyResample_clear_transformations(void* self)
{
    try {
        NiftyResample<float>& res = objectFromHandle<NiftyResample<float> >(self);
        res.clear_transformations();
        return new DataHandle;
    }
    CATCH;
}

extern "C"
void* cReg_NiftyResample_process(void* ptr)
{
    try {
        NiftyResample<float>& res = objectFromHandle<NiftyResample<float> >(ptr);
        res.process();
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftyResample_forward(const void* output_ptr, const void * const input_ptr, const void * resampler_ptr)
{
    try {
        // Get resampler
        std::shared_ptr<NiftyResample<float> > resampler_sptr;
        getObjectSptrFromHandle<NiftyResample<float> >(resampler_ptr, resampler_sptr);

        // Get input and output images
        std::shared_ptr<const ImageData> input_sptr;
        getObjectSptrFromHandle<const ImageData>(input_ptr, input_sptr);
        std::shared_ptr<ImageData> output_sptr;
        getObjectSptrFromHandle<ImageData>(output_ptr, output_sptr);

        // Forward transformation
        resampler_sptr->forward(output_sptr,input_sptr);

        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_NiftyResample_adjoint(const void* output_ptr, const void * const input_ptr, const void * resampler_ptr)
{
    try {
        // Get resampler
        std::shared_ptr<NiftyResample<float> > resampler_sptr;
        getObjectSptrFromHandle<NiftyResample<float> >(resampler_ptr, resampler_sptr);

        // Get input and output images
        std::shared_ptr<const ImageData> input_sptr;
        getObjectSptrFromHandle<const ImageData>(input_ptr, input_sptr);
        std::shared_ptr<ImageData> output_sptr;
        getObjectSptrFromHandle<ImageData>(output_ptr, output_sptr);

        // Forward transformation
        resampler_sptr->adjoint(output_sptr,input_sptr);

        return new DataHandle;
    }
    CATCH;
}
// -------------------------------------------------------------------------------- //
//      ImageWeightedMean
// -------------------------------------------------------------------------------- //
extern "C"
void* cReg_ImageWeightedMean_add_image(void* ptr, const void *obj, const float weight)
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
void* cReg_ImageWeightedMean_add_image_filename(void* ptr, const char* filename, const float weight)
{
    try {
        ImageWeightedMean<float>& im_weight = objectFromHandle<ImageWeightedMean<float> >(ptr);
        im_weight.add_image(NiftiImageData<float>(filename),weight);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_ImageWeightedMean_process(void* ptr)
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
void* cReg_Transformation_get_as_deformation_field(const void* ptr, const char* name, const void* ref)
{
    try {
        Transformation<float> *trans;

        if (strcmp(name,"AffineTransformation") == 0)
            trans = &objectFromHandle<AffineTransformation<float> >(ptr);
        else if (strcmp(name,"NiftiImageData3DDisplacement") == 0)
            trans = &objectFromHandle<NiftiImageData3DDisplacement<float> >(ptr);
        else if (strcmp(name,"NiftiImageData3DDeformation") == 0)
            trans = &objectFromHandle<NiftiImageData3DDeformation<float> >(ptr);
        else
            throw std::runtime_error("cReg_Transformation_get_as_deformation_field: type should be affine, disp or def.");

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
void* cReg_AffineTransformation_construct_from_TM(size_t ptr_TM)
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
void* cReg_AffineTransformation_construct_from_trans_and_quaternion(size_t trans_ptr, const void* quat_ptr)
{
    try {
        Quaternion<float>& quat = objectFromHandle<Quaternion<float> >(quat_ptr);
        std::array<float,3> trans;
        for (unsigned i=0; i<3; ++i)
            trans[i] = ((float*)trans_ptr)[i];
        return newObjectHandle(
                    std::make_shared<AffineTransformation<float> >(trans,quat));
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_construct_from_trans_and_euler(size_t trans_ptr, size_t euler_ptr)
{
    try {
        std::array<float,3> trans, euler;
        for (unsigned i=0; i<3; ++i) {
            trans[i] = ((float*)trans_ptr)[i];
            euler[i] = ((float*)euler_ptr)[i];
        }
        return newObjectHandle(
                    std::make_shared<AffineTransformation<float> >(trans,euler));
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_deep_copy(const void* ptr)
{
    try {
        AffineTransformation<float>& mat = objectFromHandle<AffineTransformation<float> >(ptr);
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>(mat.deep_copy()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_write(const void* ptr, const char* filename)
{
    try {
        AffineTransformation<float>& mat = objectFromHandle<AffineTransformation<float> >(ptr);
        mat.write(filename);
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_as_array(const void* ptr, size_t ptr_TM)
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
void* cReg_AffineTransformation_get_identity()
{
    try {
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>);
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_get_inverse(const void* ptr)
{
    try {
        AffineTransformation<float>& tm = objectFromHandle<AffineTransformation<float> >(ptr);
        std::shared_ptr<AffineTransformation<float> > sptr(new AffineTransformation<float>(tm.get_inverse()));
        return newObjectHandle(sptr);
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_get_Euler_angles(const void* ptr, size_t Euler)
{
    try {
        AffineTransformation<float>& tm = objectFromHandle<AffineTransformation<float> >(ptr);
        std::array<float,3> Euler_array = tm.get_Euler_angles();
        float* Euler_float = (float*)Euler;
        for (unsigned i=0; i<3; ++i)
            Euler_float[i] = Euler_array[i];
        return new DataHandle;
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_get_quaternion(const void* ptr)
{
    try {
        AffineTransformation<float>& tm = objectFromHandle<AffineTransformation<float> >(ptr);
        return newObjectHandle(
                    std::make_shared<Quaternion<float> >(tm.get_quaternion()));
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_mul(const void* mat1_ptr, const void* mat2_ptr)
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
void* cReg_AffineTransformation_equal(const void* mat1_ptr, const void* mat2_ptr)
{
    try {
        AffineTransformation<float>& mat1 = objectFromHandle<AffineTransformation<float> >(mat1_ptr);
        AffineTransformation<float>& mat2 = objectFromHandle<AffineTransformation<float> >(mat2_ptr);
        return dataHandle<int>(mat1 == mat2);
    }
    CATCH;
}
extern "C"
void* cReg_AffineTransformation_get_average(const void* handle_vector_ptr)
{
    const DataHandleVector handle_vector = objectFromHandle<const DataHandleVector>(handle_vector_ptr);
    std::vector<AffineTransformation<float> > vec;
    for (unsigned i=0; i<handle_vector.size(); ++i)
        vec.push_back(objectFromHandle<AffineTransformation<float> >(handle_vector.at(i)));

    return newObjectHandle(
                std::make_shared<AffineTransformation<float> >(AffineTransformation<float>::get_average(vec)));
}
extern "C"
void* cReg_Quaternion_construct_from_array(size_t arr)
{
    float* arr_float = (float*)arr;
    return newObjectHandle(
                std::make_shared<Quaternion<float> >(arr_float[0],arr_float[1],arr_float[2],arr_float[3]));
}

extern "C"
void* cReg_Quaternion_construct_from_AffineTransformation(const void* ptr)
{
    AffineTransformation<float>& tm = objectFromHandle<AffineTransformation<float> >(ptr);
    return newObjectHandle(
                std::make_shared<Quaternion<float> >(tm.get_quaternion()));
}

extern "C"
void* cReg_Quaternion_get_average(const void *handle_vector_ptr)
{
    const DataHandleVector handle_vector = objectFromHandle<const DataHandleVector>(handle_vector_ptr);
    std::vector<Quaternion<float> > vec;
    for (unsigned i=0; i<handle_vector.size(); ++i)
        vec.push_back(objectFromHandle<Quaternion<float> >(handle_vector.at(i)));

    return newObjectHandle(
                std::make_shared<Quaternion<float> >(Quaternion<float>::get_average(vec)));
}

extern "C"
void* cReg_Quaternion_as_array(const void* ptr, size_t arr)
{
    Quaternion<float>& quat = objectFromHandle<Quaternion<float> >(ptr);
    float* arr_float = (float*)arr;
    std::array<float,4> quat_data = quat.get_data();
    for (unsigned i=0; i<4; ++i)
        arr_float[i] = quat_data[i];
    return new DataHandle;
}