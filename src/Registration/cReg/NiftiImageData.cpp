/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
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

/*!
\file
\ingroup Registration
\brief Base class for SIRF image data.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData.h"
#include <nifti1_io.h>
#include "_reg_resampling.h"
#include <boost/filesystem.hpp>
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DTensor.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftyResample.h"
#include <iomanip>
#include <cmath>
#include <numeric>

// Remove NiftyReg's definition of isnan
#undef isnan

using namespace sirf;

static void check_folder_exists_if_not_create(const std::string &path)
{
    // If the folder doesn't exist, create it
    boost::filesystem::path path_boost(path);
    if (!boost::filesystem::exists(path_boost.parent_path())) {
        if (path_boost.parent_path().string() != "") {
            std::cout << "\n\tCreating folder: \"" << path_boost.parent_path().string() << "\"\n" << std::flush;
            boost::filesystem::create_directory(path_boost.parent_path());
        }
    }
}

template<class dataType>
NiftiImageData<dataType>::NiftiImageData(const NiftiImageData<dataType>& to_copy)
{
    *this = dynamic_cast<const ImageData&>(to_copy);
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator=(const NiftiImageData<dataType>& to_copy)
{
    *this = dynamic_cast<const ImageData&>(to_copy);
    return *this;
}

template<class dataType>
NiftiImageData<dataType>::NiftiImageData(const ImageData& to_copy)
{
    *this = to_copy;
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator=(const ImageData& to_copy)
{
    // Check for self-assignment
    if (this != &to_copy) {
        // Try to cast to NiftiImageData.
        const NiftiImageData<dataType> * const nii_ptr = dynamic_cast<const NiftiImageData<dataType> * const >(&to_copy);
        if (nii_ptr) {
            // Check the image is copyable
            if (!nii_ptr->is_initialised())
                throw std::runtime_error("Trying to copy an uninitialised image.");

            copy_nifti_image(_nifti_image,nii_ptr->_nifti_image);
            this->_data = static_cast<float*>(_nifti_image->data);
            this->_original_datatype = nii_ptr->_original_datatype;
            set_up_geom_info();
        }
        else {
            this->_nifti_image = NiftiImageData<float>::create_from_geom_info(*to_copy.get_geom_info_sptr());
            // Always float
            this->set_up_data(NIFTI_TYPE_FLOAT32);
            // Finally, copy the data
            this->copy(to_copy.begin(), this->begin(), this->end());
        }
    }
    return *this;
}

template<class dataType>
NiftiImageData<dataType>::NiftiImageData(const std::string &filename)
{
    open_nifti_image(_nifti_image,filename);
    set_up_data(_nifti_image->datatype);
}

template<class dataType>
NiftiImageData<dataType>::NiftiImageData(const nifti_image &image_nifti)
{
    copy_nifti_image(_nifti_image,std::make_shared<nifti_image>(image_nifti));
    reg_checkAndCorrectDimension(_nifti_image.get());
    set_up_data(_nifti_image->datatype);
}

template<class dataType>
std::shared_ptr<nifti_image> NiftiImageData<dataType>::create_from_geom_info(const VoxelisedGeometricalInfo3D &geom, const bool is_tensor)
{
    typedef VoxelisedGeometricalInfo3D Info;
    Info::Size            size    = geom.get_size();
    Info::Spacing         spacing = geom.get_spacing();
    Info::TransformMatrix tm      = geom.calculate_index_to_physical_point_matrix();

    int dims[8];
    dims[0] = 3;
    dims[1] = int(size[0]);
    dims[2] = int(size[1]);
    dims[3] = int(size[2]);
    dims[4] = 1;
    dims[5] = 1;
    dims[6] = 1;
    dims[7] = 1;

    // If tensor image, dims[0] and dims[5] should be 5 and 3, respectively
    if (is_tensor) {
        dims[0] = 5;
        dims[5] = 3;
    }

    nifti_image *im = nifti_make_new_nim(dims, DT_FLOAT32, 1);
    std::shared_ptr<nifti_image> _nifti_image = std::shared_ptr<nifti_image>(im, nifti_image_free);

    // Spacing
    _nifti_image->pixdim[1]=_nifti_image->dx=spacing[0];
    _nifti_image->pixdim[2]=_nifti_image->dy=spacing[1];
    _nifti_image->pixdim[3]=_nifti_image->dz=spacing[2];
    _nifti_image->pixdim[4]=0.F;
    _nifti_image->pixdim[5]=0.F;
    _nifti_image->pixdim[6]=0.F;
    _nifti_image->pixdim[7]=0.F;
    // Distances in mm
    _nifti_image->xyz_units=2;
    // Set the transformation matrix information
    _nifti_image->qform_code=1;
    AffineTransformation<float> tm_orig;
    for (unsigned i=0;i<4;++i)
        for (unsigned j=0;j<4;++j)
            tm_orig[i][j]=tm[i][j];

    AffineTransformation<float> tm_flip;
    tm_flip[0][0] = tm_flip[1][1] = -1.F; // VoxelisedGeometricalInfo3D is LPS, but nifti is RAS so flip first and second dims.
    AffineTransformation<float> tm_final = tm_flip*tm_orig;
    for (unsigned i=0;i<4;++i)
        for (unsigned j=0;j<4;++j)
            _nifti_image->qto_xyz.m[i][j]=tm_final[i][j];

    _nifti_image->qto_ijk =
            nifti_mat44_inverse(_nifti_image->qto_xyz);
    nifti_mat44_to_quatern( _nifti_image->qto_xyz,
                            &_nifti_image->quatern_b,
                            &_nifti_image->quatern_c,
                            &_nifti_image->quatern_d,
                            &_nifti_image->qoffset_x,
                            &_nifti_image->qoffset_y,
                            &_nifti_image->qoffset_z,
                            nullptr,
                            nullptr,
                            nullptr,
                            &_nifti_image->qfac );
    _nifti_image->pixdim[0]=_nifti_image->qfac;

    // Check everything is ok
    reg_checkAndCorrectDimension(_nifti_image.get());

    return _nifti_image;
}

template<class dataType>
void NiftiImageData<dataType>::construct_NiftiImageData_from_complex_im_real_component(std::shared_ptr<NiftiImageData> &out_sptr, const std::shared_ptr<const ImageData> in_sptr)
{
    // Create image from input
    out_sptr = std::make_shared<NiftiImageData<dataType> >(*in_sptr);

    auto &it_in = in_sptr->begin();
    auto &it_out = out_sptr->begin();
    for (; it_in!=in_sptr->end(); ++it_in, ++it_out)
        *it_out = (*it_in).complex_float().real();
}

template<class dataType>
void NiftiImageData<dataType>::construct_NiftiImageData_from_complex_im_imag_component(std::shared_ptr<NiftiImageData> &out_sptr, const std::shared_ptr<const ImageData> in_sptr)
{
    if (!in_sptr->is_complex())
        std::cout << "\nNiftiImageData<dataType>::construct_NiftiImageData_from_complex_im. Warning, input image is not complex. Complex component will be empty\n";

    // Create image from input
    out_sptr = std::make_shared<NiftiImageData<dataType> >(*in_sptr);

    auto &it_in = in_sptr->begin();
    auto &it_out = out_sptr->begin();
    for (; it_in!=in_sptr->end(); ++it_in, ++it_out)
        *it_out = (*it_in).complex_float().imag();
}

template<class dataType>
void NiftiImageData<dataType>::construct_NiftiImageData_from_complex_im(std::shared_ptr<NiftiImageData> &out_real_sptr, std::shared_ptr<NiftiImageData> &out_imag_sptr, const std::shared_ptr<const ImageData> in_sptr)
{
    construct_NiftiImageData_from_complex_im_real_component(out_real_sptr,in_sptr);
    construct_NiftiImageData_from_complex_im_imag_component(out_imag_sptr,in_sptr);
}

template<class dataType>
bool NiftiImageData<dataType>::operator==(const NiftiImageData<dataType> &other) const
{
    if (this == &other)
        return true;
    return this->are_equal_to_given_accuracy(*this,other,1.e-3F);
}

template<class dataType>
bool NiftiImageData<dataType>::operator!=(const NiftiImageData<dataType> &other) const
{
    return !(*this == other);
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator+=(const NiftiImageData<dataType>& rhs)
{
    maths(rhs, add);
    return *this;
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator-=(const NiftiImageData<dataType>& rhs)
{
    maths(rhs, sub);
    return *this;
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator+=(const float val)
{
    maths(val, add);
    return *this;
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator-=(const float val)
{
    maths(val, sub);
    return *this;
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator*=(const float val)
{
    maths(val, mul);
    return *this;
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator/=(const float val)
{
    maths(1.f/val, add);
    return *this;
}

template<class dataType>
float NiftiImageData<dataType>::operator()(const int index) const
{
    assert(this->is_in_bounds(index));
    return _data[index];
}

template<class dataType>
float &NiftiImageData<dataType>::operator()(const int index)
{
    assert(this->is_in_bounds(index));
    return _data[index];
}

template<class dataType>
float NiftiImageData<dataType>::operator()(const int index[7]) const
{
    assert(this->is_in_bounds(index));
    const int index_1d = this->get_1D_index(index);
    return _data[index_1d];
}

template<class dataType>
float &NiftiImageData<dataType>::operator()(const int index[7])
{
    assert(this->is_in_bounds(index));
    const int index_1d = this->get_1D_index(index);
    return _data[index_1d];
}

template<class dataType>
std::shared_ptr<const nifti_image> NiftiImageData<dataType>::get_raw_nifti_sptr() const
{
    if (!_nifti_image)
        throw std::runtime_error("Warning, nifti has not been initialised.");
    return _nifti_image;
}

template<class dataType>
std::shared_ptr<nifti_image> NiftiImageData<dataType>::get_raw_nifti_sptr()
{
    if (!_nifti_image)
        throw std::runtime_error("Warning, nifti has not been initialised.");
    return _nifti_image;
}

template<class dataType>
void NiftiImageData<dataType>::write(const std::string &filename, const int datatype) const
{
    if (!this->is_initialised())
        throw std::runtime_error("Cannot save image to file, image not initialised.");

    std::cout << "\nSaving image to file (" << filename << ")..." << std::flush;

    check_folder_exists_if_not_create(filename);

    if (_original_datatype == -1)
        throw std::runtime_error("Original datatype was not set.");

    // Create a deep copy in case we need to change datatype
    NiftiImageData copy = *this;

    // If user wants to save in a different datatype
    if (datatype != -1)
        copy.change_datatype(datatype);
    else
        copy.change_datatype(_original_datatype);

    nifti_set_filenames(copy._nifti_image.get(), filename.c_str(), 0, 0);
    nifti_image_write(copy._nifti_image.get());
    std::cout << "done.\n\n";
}

template<class dataType>
float NiftiImageData<dataType>::get_max() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_max(): Image not initialised.");

    // Get data
    return *std::max_element(_data, _data + _nifti_image->nvox);
}

template<class dataType>
float NiftiImageData<dataType>::get_min() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_min(): Image not initialised.");

    // Get data
    return *std::min_element(_data, _data + _nifti_image->nvox);
}

template<class dataType>
float NiftiImageData<dataType>::get_mean() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_min(): Image not initialised.");

    float sum = this->get_sum();
    unsigned non_nan_count = unsigned(this->get_num_voxels()) - this->get_nan_count();
    return sum / float(non_nan_count);
}

template<class dataType>
float NiftiImageData<dataType>::get_variance() const
{
    float var  = 0.f;
    float mean = this->get_mean();
    for (unsigned i=0; i<this->get_num_voxels(); ++i)
        var += std::pow(float((*this)(i)) - mean, 2.f);
    return var / float(this->get_num_voxels());
}


template<class dataType>
float NiftiImageData<dataType>::get_standard_deviation() const
{
    return sqrt(this->get_variance());
}

template<class dataType>
float NiftiImageData<dataType>::get_sum() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_sum(): Image not initialised.");

    double sum = 0;
    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        sum += double(_data[i]);
    return float(sum);
}

template<class dataType>
unsigned NiftiImageData<dataType>::get_nan_count() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_sum(): Image not initialised.");

    unsigned nan_count = 0;
    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        if (std::isnan(_data[i]))
            ++nan_count;

    return nan_count;
}

template<class dataType>
void NiftiImageData<dataType>::fill(const float v)
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::fill(): Image not initialised.");

    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        _data[i] = v;
}

template<class dataType>
float NiftiImageData<dataType>::get_norm(const NiftiImageData<dataType>& other) const
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_norm: first image is not initialised.");
    if (!other.is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_norm: second image is not initialised.");
    for (int i=0; i<8; ++i)
        if (_nifti_image->dim[i] != other._nifti_image->dim[i])
            throw std::runtime_error("NiftiImageData<dataType>::get_norm: dimensions do not match.");

    // Use double precision to minimise rounding errors
    double result(0);
    size_t num_vox = _nifti_image->nvox;
    for (size_t i=0; i<num_vox; ++i)
        // If either value is nan, skip
        if (!std::isnan(this->operator()(i)+other(i)))
            result += double(pow( this->operator()(i) - other(i), 2));

    return float(sqrt(result));
}

template<class dataType>
const int* NiftiImageData<dataType>::get_dimensions() const
{
    return _nifti_image->dim;
}

template<class dataType>
size_t NiftiImageData<dataType>::get_num_voxels() const
{
    return _nifti_image->nvox;
}

template<class dataType>
void NiftiImageData<dataType>::check_dimensions(const NiftiImageDataType image_type)
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::check_dimensions(): Image not initialised.");

    int ndim, nt, nu, intent_code, intent_p1;
    if        (image_type == _general)  { ndim=-1; nt=-1; nu=-1; intent_code = NIFTI_INTENT_NONE;   intent_p1=-1;         }
    else   if (image_type == _3D)       { ndim= 3; nt= 1; nu= 1; intent_code = NIFTI_INTENT_NONE;   intent_p1=-1;         }
    else   if (image_type == _3DTensor) { ndim= 5; nt= 1; nu= 3; intent_code = NIFTI_INTENT_VECTOR; intent_p1=-1;         }
    else   if (image_type == _3DDisp)   { ndim= 5; nt= 1; nu= 3; intent_code = NIFTI_INTENT_VECTOR; intent_p1=DISP_FIELD; }
    else /*if (image_type == _3DDef)*/  { ndim= 5; nt= 1; nu= 3; intent_code = NIFTI_INTENT_VECTOR; intent_p1=DEF_FIELD;  }

    // Check everthing is as it should be. -1 means we don't care about it
    // (e.g., NiftiImageData3D doesn't care about intent_p1, which is used by NiftyReg for Disp/Def fields)
    bool everything_ok = true;
    if (ndim         != -1 && ndim        != _nifti_image->ndim)        everything_ok = false;
    if ( nu          != -1 && nu          != _nifti_image->nu  )        everything_ok = false;
    if ( nt          != -1 && nt          != _nifti_image->nt  )        everything_ok = false;
    if ( intent_code != -1 && intent_code != _nifti_image->intent_code) everything_ok = false;
    if ( intent_p1   != -1 && intent_p1   != _nifti_image->intent_p1)   everything_ok = false;

    if (everything_ok)
        return;

    // If not, throw an error.
    std::stringstream ss;
    ss << "Trying to construct a ";
    if      (typeid(*this) == typeid(NiftiImageData3D<dataType>))             ss << "NiftiImageData3D";
    else if (typeid(*this) == typeid(NiftiImageData3DTensor<dataType>))       ss << "NiftiImageData3DTensor";
    else if (typeid(*this) == typeid(NiftiImageData3DDisplacement<dataType>)) ss << "NiftiImageData3DDisplacement";
    else if (typeid(*this) == typeid(NiftiImageData3DDeformation<dataType>))  ss << "NiftiImageData3DDeformation";
    ss << ".\n\t\tExpected params: ndim = " << ndim << ", nu = " << nu << ", nt = " << nt;
    if      (intent_code == NIFTI_INTENT_NONE)   ss << ", intent_code = None";
    else if (intent_code == NIFTI_INTENT_VECTOR) ss << ", intent_code = Vector";
    if      (intent_p1 == 0) ss << ", intent_p1 = Deformation";
    else if (intent_p1 == 1) ss << ", intent_p1 = Displacement";
    ss << "\n\t\tActual params:   ndim = " << _nifti_image->ndim << ", nu = " << _nifti_image->nu << ", nt = " << _nifti_image->nt;
    if      (_nifti_image->intent_code == NIFTI_INTENT_NONE)   ss << ", intent_code = None";
    else if (_nifti_image->intent_code == NIFTI_INTENT_VECTOR) ss << ", intent_code = Vector";
    if      (intent_p1 != -1 && _nifti_image->intent_p1 == 0)  ss << ", intent_p1 = Deformation";
    else if (intent_p1 != -1 && _nifti_image->intent_p1 == 1)  ss << ", intent_p1 = Displacement";
    //std::cout << ss.str() << "\n";
    throw std::runtime_error(ss.str());
}

template<class dataType>
void NiftiImageData<dataType>::maths(const NiftiImageData<dataType>& c, const MathsType type)
{
    if (!this->is_initialised() || !c.is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::maths_image: at least one image is not initialised.");
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*this, c, true))
        throw std::runtime_error("NiftiImageData<dataType>::maths_image: metadata do not match.");
    if (type != add && type != sub)
        throw std::runtime_error("NiftiImageData<dataType>::maths_image: only implemented for add and subtract.");

    for (int i=0; i<int(this->_nifti_image->nvox); ++i) {
        if (type == add) (*this)(i) += c(i);
        else             (*this)(i) -= c(i);
    }
}

template<class dataType>
void NiftiImageData<dataType>::maths(const float val, const MathsType type)
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::maths_image_val: image is not initialised.");
    if (type != add && type != sub && type != mul)
        throw std::runtime_error("NiftiImageData<dataType>::maths_image_val: only implemented for add, subtract and multiply.");

    for (int i=0; i<int(this->_nifti_image->nvox); ++i) {
        if      (type == add) (*this)(i) += val;
        else if (type == sub) (*this)(i) -= val;
        else                  (*this)(i) *= val;
    }
}

/// Open nifti image
template<class dataType>
void NiftiImageData<dataType>::open_nifti_image(std::shared_ptr<nifti_image> &image, const std::string &filename)
{
    // Check filename has been entered, file exists and file is nifti
    if (filename.empty())
        throw std::runtime_error("Empty filename has been supplied, cannot open nifti image.");
    if (is_nifti_file(filename.c_str()) == -1)
        throw std::runtime_error("Attempting to open a file that is not a NIFTI image.\n\tFilename: " + filename);

    // Open file
    nifti_image *im = nifti_image_read(filename.c_str(), 1);
    image = std::shared_ptr<nifti_image>(im, nifti_image_free);

    // Ensure the image has all the values correctly set
    reg_checkAndCorrectDimension(image.get());
}

/// Copy nifti image
template<class dataType>
void NiftiImageData<dataType>::copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr)
{
#ifndef NDEBUG
    std::cout << "\nPerforming hard copy of nifti image..." << std::flush;
#endif

    // Copy the info
    nifti_image *output_ptr;

    output_ptr = nifti_copy_nim_info(image_to_copy_sptr.get());
    output_image_sptr = std::shared_ptr<nifti_image>(output_ptr, nifti_image_free);

    // How much memory do we need to copy?
    size_t mem = output_image_sptr->nvox * unsigned(output_image_sptr->nbyper);

    // Allocate the memory
    output_image_sptr->data=static_cast<void *>(malloc(mem));

    // Copy!
    memcpy(output_image_sptr->data, image_to_copy_sptr->data, mem);

    // Check everything is ok
    reg_checkAndCorrectDimension(output_image_sptr.get());

#ifndef NDEBUG
    std::cout << "done.\n\n";
#endif
}

template<class dataType>
void NiftiImageData<dataType>::normalise_zero_and_one()
{
    dataType max = this->get_max();
    dataType min = this->get_min();
    // im = (im-min) / (max-min)
    for (size_t i=0; i<this->get_num_voxels(); ++i)
        (*this)(i) = ((*this)(i)-min) / (max-min);
}

template<class dataType>
void NiftiImageData<dataType>::standardise()
{
    dataType mean = this->get_mean();
    dataType std  = this->get_standard_deviation();
    for (size_t i=0; i<this->get_num_voxels(); ++i)
        (*this)(i) = ((*this)(i) - mean) / std;
}

template<class dataType>
void NiftiImageData<dataType>::change_datatype(const int datatype)
{
    if      (datatype == DT_BINARY)   change_datatype<bool>(*this);
    else if (datatype == DT_INT8)     change_datatype<signed char>(*this);
    else if (datatype == DT_INT16)    change_datatype<signed short>(*this);
    else if (datatype == DT_INT32)    change_datatype<signed int>(*this);
    else if (datatype == DT_FLOAT32)  change_datatype<float>(*this);
    else if (datatype == DT_FLOAT64)  change_datatype<double>(*this);
    else if (datatype == DT_UINT8)    change_datatype<unsigned char>(*this);
    else if (datatype == DT_UINT16)   change_datatype<unsigned short>(*this);
    else if (datatype == DT_UINT32)   change_datatype<unsigned int>(*this);
    else if (datatype == DT_INT64)    change_datatype<signed long long>(*this);
    else if (datatype == DT_UINT64)   change_datatype<unsigned long long>(*this);
    else if (datatype == DT_FLOAT128) change_datatype<long double>(*this);
    else
        throw std::runtime_error("Trying to change to bad/unsupported datatype (" + std::to_string(datatype) + " / " + nifti_datatype_to_string(datatype) + ").");
}

/// Dump header info
template<class dataType>
void NiftiImageData<dataType>::print_header() const
{
    NiftiImageData<dataType>::print_headers({this});
}

/// Dump multiple header info
template<class dataType>
void NiftiImageData<dataType>::print_headers(const std::vector<const NiftiImageData<dataType>*> &ims)
{
    dump_headers(ims);
}

template<class dataType>
void NiftiImageData<dataType>::crop(const int min_index[7], const int max_index[7])
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::crop: Image not initialised.");

    std::shared_ptr<nifti_image> im = _nifti_image;

    // If any min values are -1, set them to 0. If any max values are -1, set them to dim[i]-1
    int min_idx[7], max_idx[7];
    for (int i=0; i<7; ++i) {
        (min_index[i] > -1) ? min_idx[i] = min_index[i] : min_idx[i] = 0;
        (max_index[i] > -1) ? max_idx[i] = max_index[i] : max_idx[i] = im->dim[i+1]-1;
    }

    // Check the min. and max. indices are in bounds.
    // Check the max. is less than image dimensions.
    // Check that min <= max.
    bool bounds_ok = true;
    if (!this->is_in_bounds(min_idx))  bounds_ok = false;
    if (!this->is_in_bounds(max_idx))  bounds_ok = false;
    for (int i=0; i<7; ++i)
        if (max_idx[i] > im->dim[i+1]) bounds_ok = false;
    for (int i=0; i<7; ++i)
        if (min_idx[i] > max_idx[i]) bounds_ok = false;
    if (!bounds_ok) {
        std::stringstream ss;
        ss << "crop_image: Bounds not ok.\n";
        ss << "\tImage dims              = (";
        for (int i=1; i<8; ++i) ss << im->dim[i] << " ";
        ss << ").\n\tMinimum requested index = (";
        for (int i=0; i<7; ++i) ss << min_idx[i] << " ";
        ss << ").\n\tMaximum requested index = (";
        for (int i=0; i<7; ++i) ss << max_idx[i] << " ";
        ss << ").\n";
        throw std::runtime_error(ss.str());
    }

    // Copy the original array
    const NiftiImageData copy = *this;

    // Set the new number of voxels
    im->dim[1] = im->nx = max_idx[0] - min_idx[0] + 1;
    im->dim[2] = im->ny = max_idx[1] - min_idx[1] + 1;
    im->dim[3] = im->nz = max_idx[2] - min_idx[2] + 1;
    im->dim[4] = im->nt = max_idx[3] - min_idx[3] + 1;
    im->dim[5] = im->nu = max_idx[4] - min_idx[4] + 1;
    im->dim[6] = im->nv = max_idx[5] - min_idx[5] + 1;
    im->dim[7] = im->nw = max_idx[6] - min_idx[6] + 1;
    im->nvox = unsigned(im->nx * im->ny * im->nz * im->nt * im->nu * im->nv * im->nw);

    // Set the number of dimensions (if it's been decreased)
    if (im->dim[0] == 7 && im->nw == 1) im->dim[0] = im->ndim = 6;
    if (im->dim[0] == 6 && im->nv == 1) im->dim[0] = im->ndim = 5;
    if (im->dim[0] == 5 && im->nu == 1) im->dim[0] = im->ndim = 4;
    if (im->dim[0] == 4 && im->nt == 1) im->dim[0] = im->ndim = 3;
    if (im->dim[0] == 3 && im->nz == 1) im->dim[0] = im->ndim = 2;
    if (im->dim[0] == 2 && im->ny == 1) im->dim[0] = im->ndim = 1;
    if (im->dim[0] == 1 && im->nx == 1) im->dim[0] = im->ndim = 0;

    // Reset the data to the correct num of voxels
    free(im->data);
    im->data = static_cast<void*>(calloc(im->nvox,size_t(im->nbyper)));
    _data    = static_cast<float*>(im->data);

    // Get the data
    float *old_data = static_cast<float*>(copy.get_raw_nifti_sptr()->data);
    float *new_data = _data;

    int new_index[7], old_index[7];
    int new_1d_idx, old_1d_idx;

    // Fill the data
    for (old_index[6]=min_idx[6]; old_index[6]<=max_idx[6]; ++old_index[6]) {
        for (old_index[5]=min_idx[5]; old_index[5]<=max_idx[5]; ++old_index[5]) {
            for (old_index[4]=min_idx[4]; old_index[4]<=max_idx[4]; ++old_index[4]) {
                for (old_index[3]=min_idx[3]; old_index[3]<=max_idx[3]; ++old_index[3]) {
                    for (old_index[2]=min_idx[2]; old_index[2]<=max_idx[2]; ++old_index[2]) {
                        for (old_index[1]=min_idx[1]; old_index[1]<=max_idx[1]; ++old_index[1]) {
                            for (old_index[0]=min_idx[0]; old_index[0]<=max_idx[0]; ++old_index[0]) {

                                for (int i=0; i<7; ++i)
                                    new_index[i] = old_index[i] - min_idx[i];

                                new_1d_idx = this->get_1D_index(new_index);
                                old_1d_idx = copy.get_1D_index(old_index);
                                new_data[new_1d_idx] = old_data[old_1d_idx];
                            }
                        }
                    }
                }
            }
        }
    }

    // If the minimum has been changed, need to alter the origin.
    for (int i=0; i<3; ++i)
        _nifti_image->qto_ijk.m[i][3] -= min_idx[i];
    _nifti_image->qto_xyz =
            nifti_mat44_inverse(_nifti_image->qto_ijk);
    nifti_mat44_to_quatern( _nifti_image->qto_xyz,
                            &_nifti_image->quatern_b,
                            &_nifti_image->quatern_c,
                            &_nifti_image->quatern_d,
                            &_nifti_image->qoffset_x,
                            &_nifti_image->qoffset_y,
                            &_nifti_image->qoffset_z,
                            nullptr,
                            nullptr,
                            nullptr,
                            &_nifti_image->qfac );
    _nifti_image->pixdim[0]=_nifti_image->qfac;

    this->set_up_data(DT_FLOAT32);
}

template<class dataType>
void NiftiImageData<dataType>::pad(const int min_index[7], const int max_index[7], const dataType val)
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::crop: Image not initialised.");

    std::shared_ptr<nifti_image> im = _nifti_image;

    // If any min or max values are -ve, set them to 0
    int min_idx[7], max_idx[7];
    for (int i=0; i<7; ++i) {
        (min_index[i] > -1) ? min_idx[i] = min_index[i] : min_idx[i] = 0;
        (max_index[i] > -1) ? max_idx[i] = max_index[i] : max_idx[i] = 0;
    }

    // Keep track of the old max (min is 0's)
    int old_max_idx[7];
    for (unsigned i=0; i<7; ++i)
        old_max_idx[i] = im->dim[i+1];

    // Copy the original array
    const NiftiImageData copy = *this;

    // Set the new number of voxels
    im->dim[1] = im->nx = im->dim[1] + max_idx[0] + min_idx[0];
    im->dim[2] = im->ny = im->dim[2] + max_idx[1] + min_idx[1];
    im->dim[3] = im->nz = im->dim[3] + max_idx[2] + min_idx[2];
    im->dim[4] = im->nt = im->dim[4] + max_idx[3] + min_idx[3];
    im->dim[5] = im->nu = im->dim[5] + max_idx[4] + min_idx[4];
    im->dim[6] = im->nv = im->dim[6] + max_idx[5] + min_idx[5];
    im->dim[7] = im->nw = im->dim[7] + max_idx[6] + min_idx[6];
    im->nvox = unsigned(im->nx * im->ny * im->nz * im->nt * im->nu * im->nv * im->nw);

    // Set the number of dimensions - largest non singleton
    im->dim[0] = im->ndim = 1;
    for (unsigned i=1; i<8; ++i)
        if (im->dim[i] > 1)
            im->dim[0] = im->ndim = int(i);

    // Reset the data to the correct num of voxels
    free(im->data);
    im->data = static_cast<void*>(calloc(im->nvox,size_t(im->nbyper)));
    _data    = static_cast<float*>(im->data);

    // Get the data
    float *old_data = static_cast<float*>(copy.get_raw_nifti_sptr()->data);
    float *new_data = _data;

    // Fill with the desired value
    for (unsigned i=0; i<im->nvox; ++i)
        _data[i] = val;

    // Replace the central part with the original image
    int old_index[7], new_index[7];
    int new_1d_idx, old_1d_idx;
    for (old_index[0]=0; old_index[0]<old_max_idx[0]; ++old_index[0]) {
        for (old_index[1]=0; old_index[1]<old_max_idx[1]; ++old_index[1]) {
            for (old_index[2]=0; old_index[2]<old_max_idx[2]; ++old_index[2]) {
                for (old_index[3]=0; old_index[3]<old_max_idx[3]; ++old_index[3]) {
                    for (old_index[4]=0; old_index[4]<old_max_idx[4]; ++old_index[4]) {
                        for (old_index[5]=0; old_index[5]<old_max_idx[5]; ++old_index[5]) {
                            for (old_index[6]=0; old_index[6]<old_max_idx[6]; ++old_index[6]) {

                                for (int i=0; i<7; ++i)
                                    new_index[i] = old_index[i] + min_idx[i];

                                new_1d_idx = this->get_1D_index(new_index);
                                old_1d_idx = copy.get_1D_index(old_index);
                                new_data[new_1d_idx] = old_data[old_1d_idx];
                            }
                        }
                    }
                }
            }
        }
    }

    // If the minimum has been changed, need to alter the origin.
    for (int i=0; i<3; ++i)
        _nifti_image->qto_ijk.m[i][3] += min_idx[i];
    _nifti_image->qto_xyz =
            nifti_mat44_inverse(_nifti_image->qto_ijk);
    nifti_mat44_to_quatern( _nifti_image->qto_xyz,
                            &_nifti_image->quatern_b,
                            &_nifti_image->quatern_c,
                            &_nifti_image->quatern_d,
                            &_nifti_image->qoffset_x,
                            &_nifti_image->qoffset_y,
                            &_nifti_image->qoffset_z,
                            nullptr,
                            nullptr,
                            nullptr,
                            &_nifti_image->qfac );
    _nifti_image->pixdim[0]=_nifti_image->qfac;

    this->set_up_data(DT_FLOAT32);
}

template<class dataType>
int NiftiImageData<dataType>::get_1D_index(const int idx[7]) const
{
    // Get dims and spacing
    int *dim = _nifti_image->dim;

    // Check it's in bounds
    for (int i=0; i<7; ++i) {
        if (idx[i]<0 || idx[i]>=dim[i+1]) {
            std::stringstream ss;
            ss << "NiftiImageData<dataType>::get_1D_index: Element out of bounds.\n";
            ss << "\tRequested = ( ";
            for (int i=0;i<7;++i) ss << idx[i] << " ";
            ss << ")\n\tBounds    = ( ";
            for (int i=0;i<7;++i) ss << dim[i+1] << " ";
            ss << ")";
            throw std::runtime_error(ss.str());
        }
    }

    int idx_1d = 0;
    idx_1d += idx[0];
    idx_1d += idx[1] * dim[1];
    idx_1d += idx[2] * dim[1] * dim[2];
    idx_1d += idx[3] * dim[1] * dim[2] * dim[3];
    idx_1d += idx[4] * dim[1] * dim[2] * dim[3] * dim[4];
    idx_1d += idx[5] * dim[1] * dim[2] * dim[3] * dim[4] * dim[5];
    idx_1d += idx[6] * dim[1] * dim[2] * dim[3] * dim[4] * dim[5] * dim[6];

    return idx_1d;
}

template<class dataType>
void NiftiImageData<dataType>::set_up_data(const int original_datatype)
{
    // Save the original datatype, we'll convert it back to this just before saving
    _original_datatype = original_datatype;

    // Display a warning that data will be lost if original was e.g., double
    if (original_datatype != NIFTI_TYPE_FLOAT32) {
        if (_nifti_image->nbyper > int(sizeof(float)))
            std::cout << "\nDecreasing number of bytes per pixel, could cause loss of accuracy.\n"
                 << "Input data type was " << nifti_datatype_to_string(original_datatype)
                 << ", converting to " << nifti_datatype_to_string(NIFTI_TYPE_FLOAT32) << ".\n";

        this->change_datatype(NIFTI_TYPE_FLOAT32);
    }

    _nifti_image->nbyper = sizeof(float);
    this->_data = static_cast<float*>(_nifti_image->data);

    // Take slope and intercept into account
    if (std::abs(_nifti_image->scl_slope-1) > 1e-4f || std::abs(_nifti_image->scl_inter) > 1e-4f) {
        for (unsigned i=0; i<this->get_num_voxels(); ++i)
            _data[i] = _nifti_image->scl_slope * _data[i] + _nifti_image->scl_inter;
        _nifti_image->scl_slope = 1.f;
        _nifti_image->scl_inter = 0.f;
    }

    // Lastly, initialise the geometrical info
    set_up_geom_info();
}

template<class dataType>
bool NiftiImageData<dataType>::is_in_bounds(const int index[7]) const
{
    for (int i=0; i<7; ++i)
        if (index[i]<0 || index[0]>=_nifti_image->dim[1])
            return false;
    return true;
}

template<class dataType>
bool NiftiImageData<dataType>::is_in_bounds(const int index) const
{
    return (index>=0 && index<int(_nifti_image->nvox));
}

template<class dataType>
bool NiftiImageData<dataType>::is_same_size(const NiftiImageData &im) const
{
    for (int i=0; i<8; ++i)
        if (_nifti_image->dim[i] != im._nifti_image->dim[i])
            return false;
    return true;
}

template<typename T>
static bool do_nifti_image_metadata_elements_match(const std::string &name, const T &elem1, const T &elem2, bool verbose)
{
    if(float(fabs(float(elem1-elem2))) < 1.e-7F)
        return true;
    if (verbose)
        std::cout << "mismatch in " << name << " , (values: " <<  elem1 << " and " << elem2 << ")\n";
    return false;
}

static bool do_nifti_image_metadata_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2, bool verbose)
{
    AffineTransformation<float> e1(elem1.m), e2(elem2.m);
    if(e1 == e2)
        return true;
    if (verbose) {
        std::cout << "mismatch in " << name << "\n";
        AffineTransformation<float>::print({e1, e2});
        std::cout << "\n";
    }
    return false;
}

/// Do nifti image metadatas match?
template<class dataType>
bool NiftiImageData<dataType>::do_nifti_image_metadata_match(const NiftiImageData &im1, const NiftiImageData &im2, bool verbose)
{
#ifndef NDEBUG
    std::cout << "\nChecking if metadata of two images match..." << std::flush;
#endif

    std::shared_ptr<const nifti_image> im1_sptr = im1.get_raw_nifti_sptr();
    std::shared_ptr<const nifti_image> im2_sptr = im2.get_raw_nifti_sptr();

    bool images_match =
            do_nifti_image_metadata_elements_match("analyze75_orient",im1_sptr->analyze75_orient,im2_sptr->analyze75_orient, verbose) &&
            do_nifti_image_metadata_elements_match("byteorder",       im1_sptr->byteorder,       im2_sptr->byteorder,        verbose) &&
            do_nifti_image_metadata_elements_match("cal_max",         im1_sptr->cal_max,         im2_sptr->cal_max,          verbose) &&
            do_nifti_image_metadata_elements_match("cal_min",         im1_sptr->cal_min,         im2_sptr->cal_min,          verbose) &&
            do_nifti_image_metadata_elements_match("datatype",        im1_sptr->datatype,        im2_sptr->datatype,         verbose) &&
            //do_nifti_image_metadata_elements_match("ext_list",        im1_sptr->ext_list,        im2_sptr->ext_list,         verbose) &&
            do_nifti_image_metadata_elements_match("freq_dim",        im1_sptr->freq_dim,        im2_sptr->freq_dim,         verbose) &&
            do_nifti_image_metadata_elements_match("iname_offset",    im1_sptr->iname_offset,    im2_sptr->iname_offset,     verbose) &&
            do_nifti_image_metadata_elements_match("intent_code",     im1_sptr->intent_code,     im2_sptr->intent_code,      verbose) &&
            do_nifti_image_metadata_elements_match("intent_p1",       im1_sptr->intent_p1,       im2_sptr->intent_p1,        verbose) &&
            do_nifti_image_metadata_elements_match("intent_p2",       im1_sptr->intent_p2,       im2_sptr->intent_p2,        verbose) &&
            do_nifti_image_metadata_elements_match("intent_p3",       im1_sptr->intent_p3,       im2_sptr->intent_p3,        verbose) &&
            do_nifti_image_metadata_elements_match("nbyper",          im1_sptr->nbyper,          im2_sptr->nbyper,           verbose) &&
            do_nifti_image_metadata_elements_match("ndim",            im1_sptr->ndim,            im2_sptr->ndim,             verbose) &&
            do_nifti_image_metadata_elements_match("nifti_type",      im1_sptr->nifti_type,      im2_sptr->nifti_type,       verbose) &&
            do_nifti_image_metadata_elements_match("nt",              im1_sptr->nt,              im2_sptr->nt,               verbose) &&
            do_nifti_image_metadata_elements_match("nu",              im1_sptr->nu,              im2_sptr->nu,               verbose) &&
            do_nifti_image_metadata_elements_match("num_ext",         im1_sptr->num_ext,         im2_sptr->num_ext,          verbose) &&
            do_nifti_image_metadata_elements_match("nv",              im1_sptr->nv,              im2_sptr->nv,               verbose) &&
            do_nifti_image_metadata_elements_match("nvox",            im1_sptr->nvox,            im2_sptr->nvox,             verbose) &&
            do_nifti_image_metadata_elements_match("nw",              im1_sptr->nw,              im2_sptr->nw,               verbose) &&
            do_nifti_image_metadata_elements_match("nx",              im1_sptr->nx,              im2_sptr->nx,               verbose) &&
            do_nifti_image_metadata_elements_match("ny",              im1_sptr->ny,              im2_sptr->ny,               verbose) &&
            do_nifti_image_metadata_elements_match("nz",              im1_sptr->nz,              im2_sptr->nz,               verbose) &&
            do_nifti_image_metadata_elements_match("phase_dim",       im1_sptr->phase_dim,       im2_sptr->phase_dim,        verbose) &&
            do_nifti_image_metadata_elements_match("qfac",            im1_sptr->qfac,            im2_sptr->qfac,             verbose) &&
            do_nifti_image_metadata_elements_match("qform_code",      im1_sptr->qform_code,      im2_sptr->qform_code,       verbose) &&
            do_nifti_image_metadata_elements_match("qoffset_x",       im1_sptr->qoffset_x,       im2_sptr->qoffset_x,        verbose) &&
            do_nifti_image_metadata_elements_match("qoffset_y",       im1_sptr->qoffset_y,       im2_sptr->qoffset_y,        verbose) &&
            do_nifti_image_metadata_elements_match("qoffset_z",       im1_sptr->qoffset_z,       im2_sptr->qoffset_z,        verbose) &&
            do_nifti_image_metadata_elements_match("quatern_b",       im1_sptr->quatern_b,       im2_sptr->quatern_b,        verbose) &&
            do_nifti_image_metadata_elements_match("quatern_c",       im1_sptr->quatern_c,       im2_sptr->quatern_c,        verbose) &&
            do_nifti_image_metadata_elements_match("quatern_d",       im1_sptr->quatern_d,       im2_sptr->quatern_d,        verbose) &&
            do_nifti_image_metadata_elements_match("scl_inter",       im1_sptr->scl_inter,       im2_sptr->scl_inter,        verbose) &&
            do_nifti_image_metadata_elements_match("scl_slope",       im1_sptr->scl_slope,       im2_sptr->scl_slope,        verbose) &&
            do_nifti_image_metadata_elements_match("sform_code",      im1_sptr->sform_code,      im2_sptr->sform_code,       verbose) &&
            do_nifti_image_metadata_elements_match("slice_code",      im1_sptr->slice_code,      im2_sptr->slice_code,       verbose) &&
            do_nifti_image_metadata_elements_match("slice_dim",       im1_sptr->slice_dim,       im2_sptr->slice_dim,        verbose) &&
            do_nifti_image_metadata_elements_match("slice_duration",  im1_sptr->slice_duration,  im2_sptr->slice_duration,   verbose) &&
            do_nifti_image_metadata_elements_match("slice_end",       im1_sptr->slice_end,       im2_sptr->slice_end,        verbose) &&
            do_nifti_image_metadata_elements_match("slice_start",     im1_sptr->slice_start,     im2_sptr->slice_start,      verbose) &&
            do_nifti_image_metadata_elements_match("swapsize",        im1_sptr->swapsize,        im2_sptr->swapsize,         verbose) &&
            do_nifti_image_metadata_elements_match("time_units",      im1_sptr->time_units,      im2_sptr->time_units,       verbose) &&
            do_nifti_image_metadata_elements_match("toffset",         im1_sptr->toffset,         im2_sptr->toffset,          verbose) &&
            do_nifti_image_metadata_elements_match("xyz_units",       im1_sptr->xyz_units,       im2_sptr->xyz_units,        verbose) &&
            do_nifti_image_metadata_elements_match("qto_ijk",         im1_sptr->qto_ijk,         im2_sptr->qto_ijk,          verbose) &&
            do_nifti_image_metadata_elements_match("qto_xyz",         im1_sptr->qto_xyz,         im2_sptr->qto_xyz,          verbose) &&
            do_nifti_image_metadata_elements_match("sto_ijk",         im1_sptr->sto_ijk,         im2_sptr->sto_ijk,          verbose) &&
            do_nifti_image_metadata_elements_match("sto_xyz",         im1_sptr->sto_xyz,         im2_sptr->sto_xyz,          verbose);

    for (int i=0; i<8; i++) {
        if (!do_nifti_image_metadata_elements_match("dim["+std::to_string(i)+"]",    im1_sptr->dim[i],    im2_sptr->dim[i],    verbose)) images_match = false;
        // only check the dimensions of non singleton dimensions
        if (i<=im1_sptr->ndim)
            if (!do_nifti_image_metadata_elements_match("pixdim["+std::to_string(i)+"]", im1_sptr->pixdim[i], im2_sptr->pixdim[i], verbose)) images_match = false;
    }

#ifndef NDEBUG
    if (images_match) std::cout << "\tOK!\n";
#endif

    return images_match;
}

/// Dump info of multiple nifti images
template<class dataType>
void NiftiImageData<dataType>::dump_headers(const std::vector<const NiftiImageData<dataType>*> &ims)
{
    std::cout << "\nPrinting info for " << ims.size() << " nifti image(s):\n";
    dump_nifti_element(ims, "analyze_75_orient", &nifti_image::analyze75_orient);
    dump_nifti_element(ims, "analyze75_orient",  &nifti_image::analyze75_orient);
    dump_nifti_element(ims, "byteorder",         &nifti_image::byteorder);
    dump_nifti_element(ims, "cal_max",           &nifti_image::cal_max);
    dump_nifti_element(ims, "cal_min",           &nifti_image::cal_min);
    dump_nifti_element(ims, "datatype",          &nifti_image::datatype);
    dump_nifti_element(ims, "dt",                &nifti_image::dt);
    dump_nifti_element(ims, "du",                &nifti_image::du);
    dump_nifti_element(ims, "dv",                &nifti_image::dv);
    dump_nifti_element(ims, "dw",                &nifti_image::dw);
    dump_nifti_element(ims, "dx",                &nifti_image::dx);
    dump_nifti_element(ims, "dy",                &nifti_image::dy);
    dump_nifti_element(ims, "dz",                &nifti_image::dz);
    dump_nifti_element(ims, "ext_list",          &nifti_image::ext_list);
    dump_nifti_element(ims, "freq_dim",          &nifti_image::freq_dim);
    dump_nifti_element(ims, "iname_offset",      &nifti_image::iname_offset);
    dump_nifti_element(ims, "intent_code",       &nifti_image::intent_code);
    dump_nifti_element(ims, "intent_p1",         &nifti_image::intent_p1);
    dump_nifti_element(ims, "intent_p2",         &nifti_image::intent_p2);
    dump_nifti_element(ims, "intent_p3",         &nifti_image::intent_p3);
    dump_nifti_element(ims, "nbyper",            &nifti_image::nbyper);
    dump_nifti_element(ims, "ndim",              &nifti_image::ndim);
    dump_nifti_element(ims, "nifti_type",        &nifti_image::nifti_type);
    dump_nifti_element(ims, "num_ext",           &nifti_image::num_ext);
    dump_nifti_element(ims, "nvox",              &nifti_image::nvox);
    dump_nifti_element(ims, "nx",                &nifti_image::nx);
    dump_nifti_element(ims, "ny",                &nifti_image::ny);
    dump_nifti_element(ims, "nz",                &nifti_image::nz);
    dump_nifti_element(ims, "nt",                &nifti_image::nt);
    dump_nifti_element(ims, "nu",                &nifti_image::nu);
    dump_nifti_element(ims, "nv",                &nifti_image::nv);
    dump_nifti_element(ims, "nw",                &nifti_image::nw);
    dump_nifti_element(ims, "phase_dim",         &nifti_image::phase_dim);
    dump_nifti_element(ims, "qfac",              &nifti_image::qfac);
    dump_nifti_element(ims, "qform_code",        &nifti_image::qform_code);
    dump_nifti_element(ims, "qoffset_x",         &nifti_image::qoffset_x);
    dump_nifti_element(ims, "qoffset_y",         &nifti_image::qoffset_y);
    dump_nifti_element(ims, "qoffset_z",         &nifti_image::qoffset_z);
    dump_nifti_element(ims, "quatern_b",         &nifti_image::quatern_b);
    dump_nifti_element(ims, "quatern_c",         &nifti_image::quatern_c);
    dump_nifti_element(ims, "quatern_d",         &nifti_image::quatern_d);
    dump_nifti_element(ims, "scl_inter",         &nifti_image::scl_inter);
    dump_nifti_element(ims, "scl_slope",         &nifti_image::scl_slope);
    dump_nifti_element(ims, "sform_code",        &nifti_image::sform_code);
    dump_nifti_element(ims, "slice_code",        &nifti_image::slice_code);
    dump_nifti_element(ims, "slice_dim",         &nifti_image::slice_dim);
    dump_nifti_element(ims, "slice_duration",    &nifti_image::slice_duration);
    dump_nifti_element(ims, "slice_end",         &nifti_image::slice_end);
    dump_nifti_element(ims, "slice_start",       &nifti_image::slice_start);
    dump_nifti_element(ims, "swapsize",          &nifti_image::swapsize);
    dump_nifti_element(ims, "time_units",        &nifti_image::time_units);
    dump_nifti_element(ims, "toffset",           &nifti_image::toffset);
    dump_nifti_element(ims, "xyz_units",         &nifti_image::xyz_units);
    dump_nifti_element(ims, "dim",               &nifti_image::dim,    8);
    dump_nifti_element(ims, "pixdim",            &nifti_image::pixdim, 8);

    std::vector<std::shared_ptr<const nifti_image> > images;
    for(unsigned i=0;i<ims.size();i++)
        images.push_back(ims[i]->get_raw_nifti_sptr());

    // Print transformation matrices
    std::vector<AffineTransformation<float> > qto_ijk_vec, qto_xyz_vec, sto_ijk_vec, sto_xyz_vec;
    for(unsigned j=0; j<images.size(); j++) {
        qto_ijk_vec.push_back(images[j]->qto_ijk.m);
        qto_xyz_vec.push_back(images[j]->qto_xyz.m);
        sto_ijk_vec.push_back(images[j]->sto_ijk.m);
        sto_xyz_vec.push_back(images[j]->sto_xyz.m);
    }
    std::cout << "\t" << std::left << std::setw(19) << "qto_ijk:" << "\n";
    AffineTransformation<float>::print(qto_ijk_vec);
    std::cout << "\t" << std::left << std::setw(19) << "qto_xyz:" << "\n";
    AffineTransformation<float>::print(qto_xyz_vec);
    std::cout << "\t" << std::left << std::setw(19) << "sto_ijk:" << "\n";
    AffineTransformation<float>::print(sto_ijk_vec);
    std::cout << "\t" << std::left << std::setw(19) << "sto_xyz:" << "\n";
    AffineTransformation<float>::print(sto_xyz_vec);

    // Print original datatype
    std::string original_datatype = "orig_datatype: ";
    std::cout << "\t" << std::left << std::setw(19) << original_datatype;
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i]->get_original_datatype();

    // Print min
    std::string min_header = "min: ";
    std::cout << "\n\t" << std::left << std::setw(19) << min_header;
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i]->get_min();

    // Print max
    std::cout << "\n\t" << std::left << std::setw(19) << "max: ";
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i]->get_max();

    // Print mean
    std::cout << "\n\t" << std::left << std::setw(19) << "mean: ";
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i]->get_mean();

    // Print if image contains nans
    std::cout << "\n\t" << std::left << std::setw(19) << "contains nans?: ";
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i]->get_contains_nans();

    std::cout << "\n\n";
}

template<class dataType>
template<typename T>
void NiftiImageData<dataType>::dump_nifti_element(const std::vector<const NiftiImageData*> &ims, const std::string &name, const T &call_back)
{
    std::string header = name + ": ";
    std::cout << "\t" << std::left << std::setw(19) << header;
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i]->get_raw_nifti_sptr().get()->*call_back;
    std::cout << "\n";
}

template<class dataType>
template<typename T>
void NiftiImageData<dataType>::dump_nifti_element(const std::vector<const NiftiImageData*> &ims, const std::string &name, const T &call_back, const unsigned num_elems)
{
    for(unsigned i=0; i<num_elems; i++) {
        std::string header = name + "[" + std::to_string(i) + "]: ";
        std::cout << "\t" << std::left << std::setw(19) << header;
        for(unsigned j=0; j<ims.size(); j++)
            std::cout << std::setw(19) << (ims[j]->get_raw_nifti_sptr().get()->*call_back)[i];
        std::cout << "\n";
    }
}

template<class dataType>
void NiftiImageData<dataType>::set_voxel_spacing(const float new_spacing[3], const int interpolation_order)
{
#ifndef NDEBUG
    std::cout << "\nResampling image from voxel sizes of (" << _nifti_image->dx << ", " << _nifti_image->dy << ", " << _nifti_image->dz << ") to "
                 "(" << new_spacing[0] << ", " << new_spacing[1] << ", " << new_spacing[2] << ")\n";
#endif

    // Check image has been initialised
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::set_voxel_spacing: Image not initialised.");

    // Check all spacings are > 0
    for (int i=0; i<3; ++i)
        if (new_spacing[i] <= 0.F)
            throw std::runtime_error("NiftiImageData<dataType>::set_voxel_spacing(): New spacings must be > 0.");

    // If no changes, return
    if (std::abs(new_spacing[0]-_nifti_image->dx) < 1E-4F && std::abs(new_spacing[1]-_nifti_image->dx) < 1E-4F && std::abs(new_spacing[2]-_nifti_image->dx) < 1E-4F)
        return;

    // Check interpolation order is 0, 1 or 3.
    if (interpolation_order != 0 && interpolation_order != 1 && interpolation_order != 3)
        throw std::runtime_error("NiftiImageData<dataType>::set_voxel_spacing(): Interpolation order should be 0, 1 or 3 (NN, linear or cubic, respectively.");

    // Define the size of the new image
    int newDim[8];
    for(size_t i=0; i<8; ++i) newDim[i]=_nifti_image->dim[i];

    newDim[1]=int(ceilf(float(_nifti_image->dim[1])*_nifti_image->pixdim[1]/new_spacing[0]));
    newDim[2]=int(ceilf(float(_nifti_image->dim[2])*_nifti_image->pixdim[2]/new_spacing[1]));
    if(_nifti_image->nz>1)
        newDim[3]=int(ceilf(float(_nifti_image->dim[3])*_nifti_image->pixdim[3]/new_spacing[2]));

    // Create copy of old image
    NiftiImageData<dataType> old = *this;
    nifti_image *oldImg = old.get_raw_nifti_sptr().get();
    // Create the new image
    _nifti_image.reset(nifti_make_new_nim(newDim,_nifti_image->datatype,true),nifti_image_free);
    nifti_image *newImg = _nifti_image.get();

    newImg->pixdim[1]=newImg->dx=new_spacing[0];
    newImg->pixdim[2]=newImg->dy=new_spacing[1];
    if(oldImg->nz>1)
        newImg->pixdim[3]=newImg->dz=new_spacing[2];
    newImg->qform_code=oldImg->qform_code;
    newImg->sform_code=oldImg->sform_code;
    // Update the qform matrix
    newImg->qfac=oldImg->qfac;
    newImg->quatern_b=oldImg->quatern_b;
    newImg->quatern_c=oldImg->quatern_c;
    newImg->quatern_d=oldImg->quatern_d;
    newImg->qoffset_x=oldImg->qoffset_x+newImg->dx/2.f-oldImg->dx/2.f;
    newImg->qoffset_y=oldImg->qoffset_y+newImg->dy/2.f-oldImg->dy/2.f;
    if(oldImg->nz>1)
        newImg->qoffset_z=oldImg->qoffset_z+newImg->dz/2.f-oldImg->dz/2.f;
    else newImg->qoffset_z=oldImg->qoffset_z;
    newImg->qto_xyz=nifti_quatern_to_mat44(newImg->quatern_b,
                                           newImg->quatern_c,
                                           newImg->quatern_d,
                                           newImg->qoffset_x,
                                           newImg->qoffset_y,
                                           newImg->qoffset_z,
                                           newImg->pixdim[1],
                                           newImg->pixdim[2],
                                           newImg->pixdim[3],
                                           newImg->qfac);
    newImg->qto_ijk=nifti_mat44_inverse(newImg->qto_xyz);
    if(newImg->sform_code>0) {
        // Compute the new sform
        float scalingRatio[3];
        scalingRatio[0]= newImg->dx / oldImg->dx;
        scalingRatio[1]= newImg->dy / oldImg->dy;
        if(oldImg->nz>1)
            scalingRatio[2]= newImg->dz / oldImg->dz;
        else scalingRatio[2]=1.f;
        newImg->sto_xyz.m[0][0]=oldImg->sto_xyz.m[0][0] * scalingRatio[0];
        newImg->sto_xyz.m[1][0]=oldImg->sto_xyz.m[1][0] * scalingRatio[0];
        newImg->sto_xyz.m[2][0]=oldImg->sto_xyz.m[2][0] * scalingRatio[0];
        newImg->sto_xyz.m[3][0]=oldImg->sto_xyz.m[3][0];
        newImg->sto_xyz.m[0][1]=oldImg->sto_xyz.m[0][1] * scalingRatio[1];
        newImg->sto_xyz.m[1][1]=oldImg->sto_xyz.m[1][1] * scalingRatio[1];
        newImg->sto_xyz.m[2][1]=oldImg->sto_xyz.m[2][1] * scalingRatio[1];
        newImg->sto_xyz.m[3][1]=oldImg->sto_xyz.m[3][1];
        newImg->sto_xyz.m[0][2]=oldImg->sto_xyz.m[0][2] * scalingRatio[2];
        newImg->sto_xyz.m[1][2]=oldImg->sto_xyz.m[1][2] * scalingRatio[2];
        newImg->sto_xyz.m[2][2]=oldImg->sto_xyz.m[2][2] * scalingRatio[2];
        newImg->sto_xyz.m[3][2]=oldImg->sto_xyz.m[3][2];
        newImg->sto_xyz.m[0][3]=oldImg->sto_xyz.m[0][3]+newImg->dx/2.f-oldImg->dx/2.f;
        newImg->sto_xyz.m[1][3]=oldImg->sto_xyz.m[1][3]+newImg->dy/2.f-oldImg->dy/2.f;
        if(oldImg->nz>1)
            newImg->sto_xyz.m[2][3]=oldImg->sto_xyz.m[2][3]+newImg->dz/2.f-oldImg->dz/2.f;
        else newImg->sto_xyz.m[2][3]=oldImg->sto_xyz.m[2][3];
        newImg->sto_xyz.m[3][3]=oldImg->sto_xyz.m[3][3];
        newImg->sto_ijk=nifti_mat44_inverse(newImg->sto_xyz);
    }
    reg_checkAndCorrectDimension(newImg);
    // Create a deformation field
    nifti_image *def=nifti_copy_nim_info(newImg);
    def->dim[0]=def->ndim=5;
    def->dim[4]=def->nt=1;
    def->pixdim[4]=def->dt=1.f;
    if(newImg->nz==1)
        def->dim[5]=def->nu=2;
    else def->dim[5]=def->nu=3;
    def->pixdim[5]=def->du=1.f;
    def->dim[6]=def->nv=1;
    def->pixdim[6]=def->dv=1.f;
    def->dim[7]=def->nw=1;
    def->pixdim[7]=def->dw=1.f;
    def->nvox = size_t(def->nx * def->ny * def->nz * def->nt * def->nu);
    def->nbyper = sizeof(float);
    def->datatype = NIFTI_TYPE_FLOAT32;
    def->data = static_cast<void *>(calloc(def->nvox,size_t(def->nbyper)));
    // Fill the deformation field with an identity transformation
    reg_getDeformationFromDisplacement(def);
    // Allocate and compute the Jacobian matrices
    mat33 *jacobian = static_cast<mat33 *>(malloc(size_t(def->nx * def->ny * def->nz) * sizeof(mat33)));
    for(size_t i=0;i<size_t(def->nx*def->ny*def->nz);++i)
        reg_mat33_eye(&jacobian[i]);

    if((newImg->pixdim[1]>oldImg->pixdim[1] ||
            newImg->pixdim[2]>oldImg->pixdim[2] ||
            newImg->pixdim[3]>oldImg->pixdim[3]) && interpolation_order != 0) {
        reg_resampleImage_PSF(oldImg,
                              newImg,
                              def,
                              NULL,
                              interpolation_order,
                              0.f,
                              jacobian,
                              0);
    }
    else{
        reg_resampleImage(oldImg,
                          newImg,
                          def,
                          NULL,
                          interpolation_order,
                          0.f);
    }
    free(jacobian);
    nifti_image_free(def);

    // Store the data and update geom info
    this->_data = static_cast<float*>(_nifti_image->data);
    set_up_geom_info();
}

template<class dataType>
void NiftiImageData<dataType>::kernel_convolution(const float sigma, NREG_CONV_KERNEL_TYPE conv_type)
{
    // Check image has been initialised
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::set_voxel_spacing: Image not initialised.");

    // Warning
    std::cout << "\n\n\n\nNiftiImageData<dataType>::kernel_convolution(): Warning, I haven't tested this at all!\n\n\n\n";

    float *sigma_t=new float[_nifti_image->nt];
    for(int i=0; i<_nifti_image->nt; ++i) sigma_t[i]=sigma; //-0.7355f?
    reg_tools_kernelConvolution(_nifti_image.get(),sigma_t,conv_type);
    delete []sigma_t;
}

enum FlipOrMirror {
    Flip,
    Mirror
};

template<class dataType>
void
flip_or_mirror(const FlipOrMirror flip_or_mirror, const unsigned axis, NiftiImageData<dataType> &im)
{
    if(!im.is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::flip_or_mirror: Image not initialised.");

    // copy original
    std::unique_ptr<NiftiImageData<dataType> > original_sptr = im.clone();

    // Get dims
    int dims[7];
    for (int i=0; i<7; ++i)
        dims[i] = im.get_dimensions()[i+1];
    int old_index[7], new_index[7];

    // Loop over
    for (old_index[0]=0; old_index[0]<dims[0]; ++old_index[0]) {
        for (old_index[1]=0; old_index[1]<dims[1]; ++old_index[1]) {
            for (old_index[2]=0; old_index[2]<dims[2]; ++old_index[2]) {
                for (old_index[3]=0; old_index[3]<dims[3]; ++old_index[3]) {
                    for (old_index[4]=0; old_index[4]<dims[4]; ++old_index[4]) {
                        for (old_index[5]=0; old_index[5]<dims[5]; ++old_index[5]) {
                            for (old_index[6]=0; old_index[6]<dims[6]; ++old_index[6]) {

                                // Copy old index
                                for (unsigned i=0; i<7; ++i)
                                    new_index[i] = old_index[i];

                                // If flipping, we switch the two indices not being flipped (i.e., x=-x, y=-y for flip about z)
                                if (flip_or_mirror == Flip) {
                                    for (unsigned i=0; i<3; ++i)
                                        new_index[i] = axis == i ? old_index[i] : dims[i] - old_index[i] - 1;
                                }

                                // If mirroring, flip the axis i.e., x=-x for mirror of x
                                else {
                                    new_index[axis] = dims[axis] - old_index[axis] - 1;
                                }

                                // Copy data
                                im(new_index) = (*original_sptr)(old_index);
                            }
                        }
                    }
                }
            }
        }
    }
}

template<class dataType>
void
NiftiImageData<dataType>::
flip_along_axis(const unsigned axis)
{
    if (axis > 2)
        throw std::runtime_error("NiftiImageData<dataType>::flip_along_axis: Axis to flip should be between 0 and 2.");

    flip_or_mirror(Flip,axis,*this);
}

template<class dataType>
void
NiftiImageData<dataType>::
mirror_along_axis(const unsigned axis)
{
    if (axis > 6)
        throw std::runtime_error("NiftiImageData<dataType>::mirror_along_axis: Axis to mirror should be between 0 and 6.");

    flip_or_mirror(Mirror,axis,*this);
}

template<class dataType>
dataType
NiftiImageData<dataType>::
get_inner_product(const NiftiImageData &other) const
{
    return std::inner_product(this->begin(),this->end(),other.begin(),dataType(0));
}

template<class dataType>
bool NiftiImageData<dataType>::are_equal_to_given_accuracy(const std::shared_ptr<const NiftiImageData> &im1_sptr, const std::shared_ptr<const NiftiImageData> &im2_sptr, const float required_accuracy_compared_to_max)
{
    if(!im1_sptr->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::are_equal_to_given_accuracy: Image 1 not initialised.");
    if(!im2_sptr->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::are_equal_to_given_accuracy: Image 2 not initialised.");

    // Check the number of dimensions match
    if(im1_sptr->get_dimensions()[0] != im2_sptr->get_dimensions()[0]) {
        std::cout << "\nImage comparison: different number of dimensions (" << im1_sptr->get_dimensions()[0] << " versus " << im2_sptr->get_dimensions()[0] << ").\n";
        return false;
    }

    // Get required accuracy compared to the image maxes
    float norm;
    float epsilon = (im1_sptr->get_max()+im2_sptr->get_max())/2.F;
    epsilon *= required_accuracy_compared_to_max;

    // If metadata match, get the norm
    if (do_nifti_image_metadata_match(*im1_sptr,*im2_sptr, false))
        norm = im1_sptr->get_norm(*im2_sptr);

    // If not, we'll have to resample
    else {
        std::cout << "\nImage comparison: metadata do not match, doing resampling...\n";
        NiftyResample<float> resample;
        resample.set_interpolation_type_to_nearest_neighbour();
        resample.set_reference_image(im1_sptr);
        resample.set_floating_image(im2_sptr);
        const std::shared_ptr<const NiftiImageData<dataType> > resampled_sptr =
                std::dynamic_pointer_cast<const NiftiImageData<dataType> >(
                    resample.forward(im2_sptr));
        norm = resampled_sptr->get_norm(*im1_sptr);
    }

    norm /= float(im1_sptr->get_num_voxels());

    if (norm <= epsilon)
        return true;

    std::cout << "\nImages are not equal (norm > epsilon).\n";
    std::cout << "\tmax1                              = " << im1_sptr->get_max() << "\n";
    std::cout << "\tmax2                              = " << im2_sptr->get_max() << "\n";
    std::cout << "\tmin1                              = " << im1_sptr->get_min() << "\n";
    std::cout << "\tmin2                              = " << im2_sptr->get_min() << "\n";
    std::cout << "\tmean1                             = " << im1_sptr->get_mean() << "\n";
    std::cout << "\tmean2                             = " << im2_sptr->get_mean() << "\n";
    std::cout << "\tstandard deviation1               = " << im1_sptr->get_standard_deviation() << "\n";
    std::cout << "\tstandard deviation2               = " << im2_sptr->get_standard_deviation() << "\n";
    std::cout << "\trequired accuracy compared to max = " << required_accuracy_compared_to_max << "\n";
    std::cout << "\tepsilon                           = " << epsilon << "\n";
    std::cout << "\tnorm/num_vox                      = " << norm << "\n";
    return false;
}

template<class dataType>
bool NiftiImageData<dataType>::are_equal_to_given_accuracy(const NiftiImageData &im1, const NiftiImageData &im2, const float required_accuracy_compared_to_max)
{
    return are_equal_to_given_accuracy(im1.clone(), im2.clone(), required_accuracy_compared_to_max);
}

// ------------------------------------------------------------------------------ //
// Pure virtual methods from ImageData
// ------------------------------------------------------------------------------ //
template<class dataType>
void NiftiImageData<dataType>::dot(const DataContainer& a_x, void* ptr) const
{
    const NiftiImageData<dataType>& x = dynamic_cast<const NiftiImageData<dataType>&>(a_x);
    assert(_nifti_image->nvox == x._nifti_image->nvox);
    double s = 0.0;
    for (unsigned i=0; i<this->_nifti_image->nvox; ++i)
        s += double(_data[i] * x._data[i]);
    float* ptr_s = static_cast<float*>(ptr);
    *ptr_s = float(s);
}

template<class dataType>
void NiftiImageData<dataType>::axpby(
    const void* ptr_a, const DataContainer& a_x,
    const void* ptr_b, const DataContainer& a_y)
{
    const float a = *static_cast<const float*>(ptr_a);
    const float b = *static_cast<const float*>(ptr_b);
    const NiftiImageData<dataType>& x = dynamic_cast<const NiftiImageData<dataType>&>(a_x);
    const NiftiImageData<dataType>& y = dynamic_cast<const NiftiImageData<dataType>&>(a_y);

    // If the result hasn't been initialised, make a clone of one of them
    if (!this->is_initialised())
        *this = *x.clone();

    assert(_nifti_image->nvox == x._nifti_image->nvox);
    assert(_nifti_image->nvox == y._nifti_image->nvox);

    for (unsigned i=0; i<this->_nifti_image->nvox; ++i)
        _data[i] = a * x._data[i] + b * y._data[i];
}

template<class dataType>
float NiftiImageData<dataType>::norm() const
{
    double s = 0.0;
    for (unsigned i=0; i<this->_nifti_image->nvox; ++i)
        s += double(_data[i]*_data[i]);
    return float(sqrt(s));
}

template<class dataType>
void NiftiImageData<dataType>::multiply
    (const DataContainer& a_x, const DataContainer& a_y)
{
    const NiftiImageData<dataType>& x = dynamic_cast<const NiftiImageData<dataType>&>(a_x);
    const NiftiImageData<dataType>& y = dynamic_cast<const NiftiImageData<dataType>&>(a_y);

    // If the result hasn't been initialised, make a clone of one of them
    if (!this->is_initialised())
        *this = *x.clone();

    assert(_nifti_image->nvox == x._nifti_image->nvox);
    assert(_nifti_image->nvox == y._nifti_image->nvox);

    for (unsigned i=0; i<this->_nifti_image->nvox; ++i)
        _data[i] = x._data[i] * y._data[i];
}

template<class dataType>
void NiftiImageData<dataType>::divide
    (const DataContainer& a_x, const DataContainer& a_y)
{
    const NiftiImageData<dataType>& x = dynamic_cast<const NiftiImageData<dataType>&>(a_x);
    const NiftiImageData<dataType>& y = dynamic_cast<const NiftiImageData<dataType>&>(a_y);

    // If the result hasn't been initialised, make a clone of one of them
    if (!this->is_initialised())
        *this = *x.clone();

    assert(_nifti_image->nvox == x._nifti_image->nvox);
    assert(_nifti_image->nvox == y._nifti_image->nvox);

    if (y.get_max() < 1.e-12F)
        THROW("division by zero in NiftiImageData::divide");

    for (unsigned i=0; i<this->_nifti_image->nvox; ++i)
        _data[i] = x._data[i] / abs(y._data[i]);
}

template<class dataType>
void NiftiImageData<dataType>::set_up_geom_info()
{
#ifndef NDEBUG
    if (_nifti_image->qform_code != 1)
        std::cout << "\nWarning: NiftiImageData<dataType>::set_up_geom_info will not be accurate, as qform != 1.\n";
    // TODO: Take care of xyz_units
    if (_nifti_image->xyz_units != 2)
        std::cout << "\nWarning: NiftiImageData<dataType>::set_up_geom_info "
                     "Only implemented for xyz_units==2 (mm). "
                     "This should be easy to add (adjust spacing "
                     "and offset by 10^).\n";
#endif

    // Number of voxels
    VoxelisedGeometricalInfo3D::Size size;
    for (unsigned i=0; i<3; ++i)
        size[i] = unsigned(_nifti_image->dim[i+1]);

    // Voxel spacing
    VoxelisedGeometricalInfo3D::Spacing spacing;
    for (unsigned i=0; i<3; ++i)
        spacing[i] = _nifti_image->pixdim[i+1];

    // VoxelisedGeometricalInfo3D is LPS, but nifti is RAS so flip first and second dims.
    AffineTransformation<float> tm_orig(_nifti_image->qto_xyz.m);
    AffineTransformation<float> tm_flip;
    tm_flip[0][0] = tm_flip[1][1] = -1.F;
    AffineTransformation<float> tm_final = tm_flip*tm_orig;

    // Offset
    VoxelisedGeometricalInfo3D::Offset offset;
    for (unsigned i=0; i<3; ++i)
        offset[i] = tm_final[i][3];

    // Transformation matrix
    VoxelisedGeometricalInfo3D::DirectionMatrix direction;
    for (unsigned i=0; i<3; ++i)
        for (unsigned j=0; j<3; ++j)
            direction[i][j] = tm_final[i][j] / spacing[j];

    // Initialise the geom info shared pointer
    this->set_geom_info(std::make_shared<VoxelisedGeometricalInfo3D>(
                            VoxelisedGeometricalInfo3D(offset,spacing,size,direction)));
}

namespace sirf {
template class NiftiImageData<float>;
}
