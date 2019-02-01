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

#include "sirf/cReg/NiftiImageData.h"
#include <nifti1_io.h>
#include <_reg_tools.h>
#include <boost/filesystem.hpp>
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/NiftiImageData3DTensor.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"
#include "sirf/cReg/AffineTransformation.h"
#include "sirf/cReg/NiftyResample.h"
#include <iomanip>
#include <cmath>


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
    copy_nifti_image(_nifti_image,to_copy._nifti_image);
    set_up_data(to_copy._original_datatype);
}

template<class dataType>
NiftiImageData<dataType>& NiftiImageData<dataType>::operator=(const NiftiImageData<dataType>& to_copy)
{
    // Check for self-assignment
    if (this != &to_copy) {
        // Check the image is copyable
        if (!to_copy.is_initialised())
            throw std::runtime_error("Trying to copy an uninitialised image.");
        // Copy
        copy_nifti_image(_nifti_image,to_copy._nifti_image);
        set_up_data(to_copy._original_datatype);
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
NiftiImageData<dataType>::NiftiImageData(const dataType * const data, const VoxelisedGeometricalInfo3D &geom)
{
    this->_nifti_image = create_from_geom_info(geom);

    // Always float
    this->set_up_data(NIFTI_TYPE_FLOAT32);

    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        this->_data[i] = dataType(data[i]);
}

template<class dataType>
std::shared_ptr<nifti_image> NiftiImageData<dataType>::create_from_geom_info(const VoxelisedGeometricalInfo3D &geom)
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
    tm_flip[0][0] = tm_flip[1][1] = -1.F;
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
NiftiImageData<dataType> NiftiImageData<dataType>::operator+(const NiftiImageData<dataType>& c) const
{
    return maths(c, add);
}

template<class dataType>
NiftiImageData<dataType> NiftiImageData<dataType>::operator-(const NiftiImageData<dataType>& c) const
{
    return maths(c, sub);
}

template<class dataType>
NiftiImageData<dataType> NiftiImageData<dataType>::operator+(const float& val) const
{
    return maths(val,add);
}

template<class dataType>
NiftiImageData<dataType> NiftiImageData<dataType>::operator-(const float& val) const
{
    return maths(val,sub);
}

template<class dataType>
NiftiImageData<dataType> NiftiImageData<dataType>::operator*(const float& val) const
{
    return maths(val,mul);
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

    float sum = 0.F;
    int nan_count = 0;
    for (int i=0; i<int(_nifti_image->nvox); ++i)
        if (!std::isnan(_data[i])) {
            sum += _data[i];
            ++nan_count;
        }

    // Get data
    return sum / float(nan_count);
}

template<class dataType>
float NiftiImageData<dataType>::get_sum() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::get_sum(): Image not initialised.");

    float sum = 0.F;
    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        sum += float(_data[i]);
    return sum;
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
NiftiImageData<dataType> NiftiImageData<dataType>::maths(const NiftiImageData<dataType>& c, const MathsType type) const
{
    if (!this->is_initialised() || !c.is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::maths_image: at least one image is not initialised.");
    if (!NiftiImageData<dataType>::do_nifti_image_metadata_match(*this, c, true))
        throw std::runtime_error("NiftiImageData<dataType>::maths_image: metadata do not match.");
    if (type != add && type != sub)
        throw std::runtime_error("NiftiImageData<dataType>::maths_image: only implemented for add and subtract.");

    NiftiImageData<dataType> res = *this;

    for (int i=0; i<int(this->_nifti_image->nvox); ++i) {
        if (type == add) res(i) += c(i);
        else             res(i) -= c(i);
    }

    return res;
}

template<class dataType>
NiftiImageData<dataType> NiftiImageData<dataType>::maths(const float val, const MathsType type) const
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::maths_image_val: image is not initialised.");
    if (type != add && type != sub && type != mul)
        throw std::runtime_error("NiftiImageData<dataType>::maths_image_val: only implemented for add, subtract and multiply.");

    NiftiImageData res = *this;
    for (int i=0; i<int(this->_nifti_image->nvox); ++i) {
        if      (type == add) res(i) += val;
        else if (type == sub) res(i) -= val;
        else                  res(i) *= val;
    }
    return res;
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
    // Check the max. is less than the min.
    bool bounds_ok = true;
    if (!this->is_in_bounds(min_idx))  bounds_ok = false;
    if (!this->is_in_bounds(max_idx))  bounds_ok = false;
    for (int i=0; i<7; ++i)
        if (max_idx[i] > im->dim[i+1]) bounds_ok = false;
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

    // TODO display a warning that data will be lost if original was e.g., double
    if (original_datatype != NIFTI_TYPE_FLOAT32) {
        if (_nifti_image->nbyper > int(sizeof(float)))
            std::cout << "\nDecreasing number of bytes per pixel, could cause loss of accuracy.\n"
                 << "Input data type was " << nifti_datatype_to_string(original_datatype)
                 << ", converting to " << nifti_datatype_to_string(NIFTI_TYPE_FLOAT32) << ".\n";

        this->change_datatype(NIFTI_TYPE_FLOAT32);
    }

    _nifti_image->nbyper = sizeof(float);
    this->_data = static_cast<float*>(_nifti_image->data);

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
    if(float(fabs(elem1-elem2)) < 1.e-7F)
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
            do_nifti_image_metadata_elements_match("du",              im1_sptr->du,              im2_sptr->du,               verbose) &&
            do_nifti_image_metadata_elements_match("dv",              im1_sptr->dv,              im2_sptr->dv,               verbose) &&
            do_nifti_image_metadata_elements_match("dw",              im1_sptr->dw,              im2_sptr->dw,               verbose) &&
            do_nifti_image_metadata_elements_match("dx",              im1_sptr->dx,              im2_sptr->dx,               verbose) &&
            do_nifti_image_metadata_elements_match("dy",              im1_sptr->dy,              im2_sptr->dy,               verbose) &&
            do_nifti_image_metadata_elements_match("dz",              im1_sptr->dz,              im2_sptr->dz,               verbose) &&
            do_nifti_image_metadata_elements_match("ext_list",        im1_sptr->ext_list,        im2_sptr->ext_list,         verbose) &&
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
    std::vector<AffineTransformation<dataType> > qto_ijk_vec, qto_xyz_vec, sto_ijk_vec, sto_xyz_vec;
    for(unsigned j=0; j<images.size(); j++) {
        qto_ijk_vec.push_back(images[j]->qto_ijk.m);
        qto_xyz_vec.push_back(images[j]->qto_xyz.m);
        sto_ijk_vec.push_back(images[j]->sto_ijk.m);
        sto_xyz_vec.push_back(images[j]->sto_xyz.m);
    }
    std::cout << "\t" << std::left << std::setw(19) << "qto_ijk:" << "\n";
    AffineTransformation<dataType>::print(qto_ijk_vec);
    std::cout << "\t" << std::left << std::setw(19) << "qto_xyz:" << "\n";
    AffineTransformation<dataType>::print(qto_xyz_vec);
    std::cout << "\t" << std::left << std::setw(19) << "sto_ijk:" << "\n";
    AffineTransformation<dataType>::print(sto_ijk_vec);
    std::cout << "\t" << std::left << std::setw(19) << "sto_xyz:" << "\n";
    AffineTransformation<dataType>::print(sto_xyz_vec);

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
bool NiftiImageData<dataType>::are_equal_to_given_accuracy(const NiftiImageData &im1, const NiftiImageData &im2, const float required_accuracy_compared_to_max)
{
    if(!im1.is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::are_equal_to_given_accuracy: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw std::runtime_error("NiftiImageData<dataType>::are_equal_to_given_accuracy: Image 2 not initialised.");

    // Check the number of dimensions match
    if(im1.get_dimensions()[0] != im2.get_dimensions()[0])
        return false;

    // Get required accuracy compared to the image maxes
    float norm;
    float epsilon = (im1.get_max()+im2.get_max())/2.F;
    epsilon *= required_accuracy_compared_to_max;

    // If metadata match, get the norm
    if (do_nifti_image_metadata_match(im1,im2, false))
        norm = im1.get_norm(im2);

    // If not, we'll have to resample
    else {
        std::cout << "\nImage comparison: metadata do not match, doing resampling...\n";
        NiftyResample<float> resample;
        resample.set_interpolation_type_to_nearest_neighbour();
        resample.set_reference_image(im1.clone());
        resample.set_floating_image(im2.clone());
        resample.process();
        norm = resample.get_output_sptr()->get_norm(im1);
    }

    if (norm < epsilon)
        return true;

    std::cout << "\nImages are not equal (norm > epsilon).\n";
    std::cout << "\tmax1                              = " << im1.get_max() << "\n";
    std::cout << "\tmax2                              = " << im2.get_max() << "\n";
    std::cout << "\tmin1                              = " << im1.get_min() << "\n";
    std::cout << "\tmin2                              = " << im2.get_min() << "\n";
    std::cout << "\trequired accuracy compared to max = " << required_accuracy_compared_to_max << "\n";
    std::cout << "\tepsilon                           = " << epsilon << "\n";
    std::cout << "\tnorm                              = " << norm << "\n";
    return false;
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
#endif

    // Number of voxels
    VoxelisedGeometricalInfo3D::Size size;
    for (int i=0; i<3; ++i)
        size[i] = unsigned(_nifti_image->dim[i+1]);

    // Voxel spacing
    VoxelisedGeometricalInfo3D::Spacing spacing;
    for (int i=0; i<3; ++i)
        spacing[i] = _nifti_image->pixdim[i+1];

    // Offset
    VoxelisedGeometricalInfo3D::Offset offset;
    for (int i=0; i<3; ++i)
        offset[i] = _nifti_image->qto_xyz.m[i][3];

    // Transformation matrix
    VoxelisedGeometricalInfo3D::DirectionMatrix direction;
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            direction[i][j] = _nifti_image->qto_xyz.m[i][j];

    // Initialise the geom info shared pointer
    _geom_info_sptr = std::make_shared<VoxelisedGeometricalInfo3D>(
                VoxelisedGeometricalInfo3D(offset,spacing,size,direction));
}

namespace sirf {
template class NiftiImageData<float>;
}
