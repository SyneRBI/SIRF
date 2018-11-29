/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

#include "NiftiImageData.h"
#include "SIRFRegMisc.h"
#include <nifti1_io.h>
#include <_reg_tools.h>
#include "stir_data_containers.h"
#include "NiftiImageData3D.h"
#include "NiftiImageData3DTensor.h"
#include "NiftiImageData3DDeformation.h"
#include "NiftiImageData3DDisplacement.h"
#include "SIRFRegAffineTransformation.h"

using namespace sirf;

NiftiImageData::NiftiImageData(const NiftiImageData& to_copy)
{
    copy_nifti_image(_nifti_image,to_copy._nifti_image);
    set_up_data(to_copy._original_datatype);
}

NiftiImageData& NiftiImageData::operator=(const NiftiImageData& to_copy)
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

NiftiImageData::NiftiImageData(const std::string &filename)
{
    open_nifti_image(_nifti_image,filename);
    set_up_data(_nifti_image->datatype);
}

NiftiImageData::NiftiImageData(const nifti_image &image_nifti)
{
    copy_nifti_image(_nifti_image,std::make_shared<nifti_image>(image_nifti));
    reg_checkAndCorrectDimension(_nifti_image.get());
    set_up_data(_nifti_image->datatype);
}

NiftiImageData::NiftiImageData(const std::shared_ptr<nifti_image> image_nifti)
{
    copy_nifti_image(_nifti_image,image_nifti);
    reg_checkAndCorrectDimension(_nifti_image.get());
    set_up_data(_nifti_image->datatype);
}

bool NiftiImageData::operator==(const NiftiImageData &other) const
{
    if (this == &other)
        return true;
    return this->are_equal_to_given_accuracy(other,1.e-3F);
}

/// Equality operator
bool NiftiImageData::operator!=(const NiftiImageData &other) const
{
    return !(*this == other);
}

NiftiImageData NiftiImageData::operator+(const NiftiImageData& c) const
{
    return maths(c, add);
}

NiftiImageData NiftiImageData::operator-(const NiftiImageData& c) const
{
    return maths(c, sub);
}

NiftiImageData NiftiImageData::operator+(const float& val) const
{
    return maths(val,add);
}

NiftiImageData NiftiImageData::operator-(const float& val) const
{
    return maths(val,sub);
}

NiftiImageData NiftiImageData::operator*(const float& val) const
{
    return maths(val,mul);
}

float NiftiImageData::operator()(const int index) const
{
    assert(this->is_in_bounds(index));
    return _data[index];
}

float &NiftiImageData::operator()(const int index)
{
    assert(this->is_in_bounds(index));
    return _data[index];
}

float NiftiImageData::operator()(const int index[7]) const
{
    assert(this->is_in_bounds(index));
    const int index_1d = this->get_1D_index(index);
    return _data[index_1d];
}

float &NiftiImageData::operator()(const int index[7])
{
    assert(this->is_in_bounds(index));
    const int index_1d = this->get_1D_index(index);
    return _data[index_1d];
}

std::shared_ptr<const nifti_image> NiftiImageData::get_raw_nifti_sptr() const
{
    if (!_nifti_image)
        throw std::runtime_error("Warning, nifti has not been initialised.");
    return _nifti_image;
}

std::shared_ptr<nifti_image> NiftiImageData::get_raw_nifti_sptr()
{
    if (!_nifti_image)
        throw std::runtime_error("Warning, nifti has not been initialised.");
    return _nifti_image;
}

void NiftiImageData::save_to_file(const std::string &filename, const int datatype) const
{
    if (!this->is_initialised())
        throw std::runtime_error("Cannot save image to file.");

    std::cout << "\nSaving image to file (" << filename << ")..." << std::flush;

    sirf::check_folder_exists(filename);

    if (_original_datatype == -1)
        throw std::runtime_error("Original datatype was not set.");

    // Create a deep copy in case we need to change datatype
    NiftiImageData copy = this->deep_copy();

    // If user wants to save in a different datatype
    if (datatype != -1)
        copy.change_datatype(datatype);
    else
        copy.change_datatype(_original_datatype);

    nifti_set_filenames(copy._nifti_image.get(), filename.c_str(), 0, 0);
    nifti_image_write(copy._nifti_image.get());
    std::cout << "done.\n\n";
}

float NiftiImageData::get_max() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData::get_max(): Image not initialised.");

    // Get data
    return *std::max_element(_data, _data + _nifti_image->nvox);
}

float NiftiImageData::get_min() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData::get_min(): Image not initialised.");

    // Get data
    return *std::min_element(_data, _data + _nifti_image->nvox);
}


float NiftiImageData::get_mean() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData::get_min(): Image not initialised.");

    float sum = 0.F;
    int nan_count = 0;
    for (int i=0; i<int(_nifti_image->nvox); ++i)
        if (!isnan(_data[i])) {
            sum += _data[i];
            ++nan_count;
        }

    // Get data
    return sum / float(nan_count);
}

float NiftiImageData::get_sum() const
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData::get_sum(): Image not initialised.");

    float sum = 0.F;
    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        sum += float(_data[i]);
    return sum;
}

void NiftiImageData::fill(const float v)
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData::fill(): Image not initialised.");

    for (unsigned i=0; i<_nifti_image->nvox; ++i)
        _data[i] = v;
}

float NiftiImageData::get_norm(const NiftiImageData& other) const
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData::get_norm: first image is not initialised.");
    if (!other.is_initialised())
        throw std::runtime_error("NiftiImageData::get_norm: second image is not initialised.");
    for (int i=0; i<8; ++i)
        if (_nifti_image->dim[i] != other._nifti_image->dim[i])
            throw std::runtime_error("NiftiImageData::get_norm: dimensions do not match.");

    // Use double precision to minimise rounding errors
    double result(0);
    for (int i=0; i<int(_nifti_image->nvox); ++i) {
        const float &val1 = this->operator()(i);
        const float &val2 = other(i);
        // If either value is nan, skip
        if (!isnan(val1+val2))
            result += double(pow( this->operator()(i) - other(i), 2));
    }
    return float(sqrt(result));
}

NiftiImageData NiftiImageData::deep_copy() const
{
    NiftiImageData copy;
    copy = *this;
    return copy;
}

const int* NiftiImageData::get_dimensions() const
{
    return _nifti_image->dim;
}

void NiftiImageData::check_dimensions(const NiftiImageDataType image_type)
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData::check_dimensions(): Image not initialised.");

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
    if      (typeid(*this) == typeid(NiftiImageData3D))             ss << "NiftiImageData3D";
    else if (typeid(*this) == typeid(NiftiImageData3DTensor))       ss << "NiftiImageData3DTensor";
    else if (typeid(*this) == typeid(NiftiImageData3DDisplacement)) ss << "NiftiImageData3DDisplacement";
    else if (typeid(*this) == typeid(NiftiImageData3DDeformation))  ss << "NiftiImageData3DDeformation";
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

NiftiImageData NiftiImageData::maths(const NiftiImageData& c, const MathsType type) const
{
    if (!this->is_initialised() || !c.is_initialised())
        throw std::runtime_error("NiftiImageData::maths_image_image: at least one image is not initialised.");
    if (!NiftiImageData::do_nifti_image_metadata_match(*this, c))
        throw std::runtime_error("NiftiImageData::maths_image_image: metadata do not match.");
    if (type != add && type != sub)
        throw std::runtime_error("NiftiImageData::maths_image_image: only implemented for add and subtract.");

    NiftiImageData res = this->deep_copy();

    for (int i=0; i<int(this->_nifti_image->nvox); ++i) {
        if (type == add) res(i) += c(i);
        else             res(i) -= c(i);
    }

    return res;
}

NiftiImageData NiftiImageData::maths(const float val, const MathsType type) const
{
    if (!this->is_initialised())
        throw std::runtime_error("NiftiImageData::maths_image_val: image is not initialised.");
    if (type != add && type != sub && type != mul)
        throw std::runtime_error("NiftiImageData::maths_image_val: only implemented for add, subtract and multiply.");

    NiftiImageData res = this->deep_copy();
    for (int i=0; i<int(this->_nifti_image->nvox); ++i) {
        if      (type == add) res(i) += val;
        else if (type == sub) res(i) -= val;
        else                  res(i) *= val;
    }
    return res;
}

/// Open nifti image
void NiftiImageData::open_nifti_image(std::shared_ptr<nifti_image> &image, const boost::filesystem::path &filename)
{
    // Check filename has been entered, file exists and file is nifti
    if (filename.empty())
        throw std::runtime_error("Empty filename has been supplied, cannot open nifti image.");
    if (!boost::filesystem::exists(filename))
        throw std::runtime_error("Cannot find the file: " + filename.string() + ".");
    if (is_nifti_file(filename.c_str()) == -1)
        throw std::runtime_error("Attempting to open a file that is not a NIFTI image.\n\tFilename: " + filename.string());

    // Open file
    nifti_image *im = nifti_image_read(filename.c_str(), 1);
    image = std::shared_ptr<nifti_image>(im, nifti_image_free);

    // Ensure the image has all the values correctly set
    reg_checkAndCorrectDimension(image.get());
}

/// Save nifti image
void NiftiImageData::save_nifti_image(NiftiImageData &image, const std::string &filename)
{
    if (!image.is_initialised())
        throw std::runtime_error("Cannot save image to file.");

    std::cout << "\nSaving image to file (" << filename << ")..." << std::flush;

    boost::filesystem::path filename_boost(filename);

    // If the folder doesn't exist, create it
    if (!boost::filesystem::exists(filename_boost.parent_path())) {
        if (filename_boost.parent_path().string() != "") {
            std::cout << "\n\tCreating folder: \"" << filename_boost.parent_path().string() << "\"\n" << std::flush;
            boost::filesystem::create_directory(filename_boost.parent_path());
        }
    }

    nifti_set_filenames(image.get_raw_nifti_sptr().get(), filename.c_str(), 0, 0);
    nifti_image_write(image.get_raw_nifti_sptr().get());
    std::cout << "done.\n\n";
}

/// Copy nifti image
void NiftiImageData::copy_nifti_image(std::shared_ptr<nifti_image> &output_image_sptr, const std::shared_ptr<nifti_image> &image_to_copy_sptr)
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

void NiftiImageData::change_datatype(const int datatype)
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

template<typename newType>
void NiftiImageData::change_datatype(NiftiImageData &im)
{
    if (im.get_raw_nifti_sptr()->datatype == DT_BINARY)   return change_datatype<newType,bool>              (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT8)     return change_datatype<newType,signed char>       (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT16)    return change_datatype<newType,signed short>      (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT32)    return change_datatype<newType,signed int>        (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT32)  return change_datatype<newType,float>             (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT64)  return change_datatype<newType,double>            (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT8)    return change_datatype<newType,unsigned char>     (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT16)   return change_datatype<newType,unsigned short>    (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT32)   return change_datatype<newType,unsigned int>      (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_INT64)    return change_datatype<newType,signed long long>  (im);
    if (im.get_raw_nifti_sptr()->datatype == DT_UINT64)   return change_datatype<newType,unsigned long long>(im);
    if (im.get_raw_nifti_sptr()->datatype == DT_FLOAT128) return change_datatype<newType,long double>       (im);

    std::stringstream ss;
    ss << "NiftImage::get_max not implemented for your data type: ";
    ss << nifti_datatype_string(im.get_raw_nifti_sptr()->datatype);
    ss << " (bytes per voxel: ";
    ss << im.get_raw_nifti_sptr()->nbyper << ").";
    throw std::runtime_error(ss.str());
}
template void NiftiImageData::change_datatype<bool>              (NiftiImageData &im);
template void NiftiImageData::change_datatype<signed char>       (NiftiImageData &im);
template void NiftiImageData::change_datatype<signed short>      (NiftiImageData &im);
template void NiftiImageData::change_datatype<signed int>        (NiftiImageData &im);
template void NiftiImageData::change_datatype<float>             (NiftiImageData &im);
template void NiftiImageData::change_datatype<double>            (NiftiImageData &im);
template void NiftiImageData::change_datatype<unsigned char>     (NiftiImageData &im);
template void NiftiImageData::change_datatype<unsigned short>    (NiftiImageData &im);
template void NiftiImageData::change_datatype<unsigned int>      (NiftiImageData &im);
template void NiftiImageData::change_datatype<signed long long>  (NiftiImageData &im);
template void NiftiImageData::change_datatype<unsigned long long>(NiftiImageData &im);
template void NiftiImageData::change_datatype<long double>       (NiftiImageData &im);

/// Dump header info
void NiftiImageData::print_header() const
{
    NiftiImageData::print_headers({*this});
}

/// Dump multiple header info
void NiftiImageData::print_headers(const std::vector<NiftiImageData> &ims)
{
    dump_headers(ims);
}

void NiftiImageData::crop(const int min_index[7], const int max_index[7])
{
    if(!this->is_initialised())
        throw std::runtime_error("NiftiImageData::crop: Image not initialised.");

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
    const NiftiImageData copy = this->deep_copy();

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

int NiftiImageData::get_1D_index(const int idx[7]) const
{
    // Get dims and spacing
    int *dim = _nifti_image->dim;

    // Check it's in bounds
    for (int i=0; i<7; ++i) {
        if (idx[i]<0 || idx[i]>=dim[i+1]) {
            std::stringstream ss;
            ss << "NiftiImageData::get_1D_index: Element out of bounds.\n";
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

void NiftiImageData::set_up_data(const int original_datatype)
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
}

bool NiftiImageData::is_in_bounds(const int index[7]) const
{
    for (int i=0; i<7; ++i)
        if (index[i]<0 || index[0]>=_nifti_image->dim[1])
            return false;
    return true;
}

bool NiftiImageData::is_in_bounds(const int index) const
{
    return (index>=0 && index<int(_nifti_image->nvox));
}

bool NiftiImageData::is_same_size(const NiftiImageData &im) const
{
    for (int i=0; i<8; ++i)
        if (_nifti_image->dim[i] != im._nifti_image->dim[i])
            return false;
    return true;
}

/// Do nifti image metadatas match?
bool NiftiImageData::do_nifti_image_metadata_match(const NiftiImageData &im1, const NiftiImageData &im2)
{
#ifndef NDEBUG
    std::cout << "\nChecking if metadata of two images match..." << std::flush;
#endif

    std::shared_ptr<const nifti_image> im1_sptr = im1.get_raw_nifti_sptr();
    std::shared_ptr<const nifti_image> im2_sptr = im2.get_raw_nifti_sptr();

    bool images_match =
            do_nifti_image_metadata_elements_match("analyze75_orient",im1_sptr->analyze75_orient,im2_sptr->analyze75_orient) &&
            do_nifti_image_metadata_elements_match("byteorder",       im1_sptr->byteorder,       im2_sptr->byteorder       ) &&
            do_nifti_image_metadata_elements_match("cal_max",         im1_sptr->cal_max,         im2_sptr->cal_max         ) &&
            do_nifti_image_metadata_elements_match("cal_min",         im1_sptr->cal_min,         im2_sptr->cal_min         ) &&
            do_nifti_image_metadata_elements_match("datatype",        im1_sptr->datatype,        im2_sptr->datatype        ) &&
            do_nifti_image_metadata_elements_match("du",              im1_sptr->du,              im2_sptr->du              ) &&
            do_nifti_image_metadata_elements_match("dv",              im1_sptr->dv,              im2_sptr->dv              ) &&
            do_nifti_image_metadata_elements_match("dw",              im1_sptr->dw,              im2_sptr->dw              ) &&
            do_nifti_image_metadata_elements_match("dx",              im1_sptr->dx,              im2_sptr->dx              ) &&
            do_nifti_image_metadata_elements_match("dy",              im1_sptr->dy,              im2_sptr->dy              ) &&
            do_nifti_image_metadata_elements_match("dz",              im1_sptr->dz,              im2_sptr->dz              ) &&
            do_nifti_image_metadata_elements_match("ext_list",        im1_sptr->ext_list,        im2_sptr->ext_list        ) &&
            do_nifti_image_metadata_elements_match("freq_dim",        im1_sptr->freq_dim,        im2_sptr->freq_dim        ) &&
            do_nifti_image_metadata_elements_match("iname_offset",    im1_sptr->iname_offset,    im2_sptr->iname_offset    ) &&
            do_nifti_image_metadata_elements_match("intent_code",     im1_sptr->intent_code,     im2_sptr->intent_code     ) &&
            do_nifti_image_metadata_elements_match("intent_p1",       im1_sptr->intent_p1,       im2_sptr->intent_p1       ) &&
            do_nifti_image_metadata_elements_match("intent_p2",       im1_sptr->intent_p2,       im2_sptr->intent_p2       ) &&
            do_nifti_image_metadata_elements_match("intent_p3",       im1_sptr->intent_p3,       im2_sptr->intent_p3       ) &&
            do_nifti_image_metadata_elements_match("nbyper",          im1_sptr->nbyper,          im2_sptr->nbyper          ) &&
            do_nifti_image_metadata_elements_match("ndim",            im1_sptr->ndim,            im2_sptr->ndim            ) &&
            do_nifti_image_metadata_elements_match("nifti_type",      im1_sptr->nifti_type,      im2_sptr->nifti_type      ) &&
            do_nifti_image_metadata_elements_match("nt",              im1_sptr->nt,              im2_sptr->nt              ) &&
            do_nifti_image_metadata_elements_match("nu",              im1_sptr->nu,              im2_sptr->nu              ) &&
            do_nifti_image_metadata_elements_match("num_ext",         im1_sptr->num_ext,         im2_sptr->num_ext         ) &&
            do_nifti_image_metadata_elements_match("nv",              im1_sptr->nv,              im2_sptr->nv              ) &&
            do_nifti_image_metadata_elements_match("nvox",            im1_sptr->nvox,            im2_sptr->nvox            ) &&
            do_nifti_image_metadata_elements_match("nw",              im1_sptr->nw,              im2_sptr->nw              ) &&
            do_nifti_image_metadata_elements_match("nx",              im1_sptr->nx,              im2_sptr->nx              ) &&
            do_nifti_image_metadata_elements_match("ny",              im1_sptr->ny,              im2_sptr->ny              ) &&
            do_nifti_image_metadata_elements_match("nz",              im1_sptr->nz,              im2_sptr->nz              ) &&
            do_nifti_image_metadata_elements_match("phase_dim",       im1_sptr->phase_dim,       im2_sptr->phase_dim       ) &&
            do_nifti_image_metadata_elements_match("qfac",            im1_sptr->qfac,            im2_sptr->qfac            ) &&
            do_nifti_image_metadata_elements_match("qform_code",      im1_sptr->qform_code,      im2_sptr->qform_code      ) &&
            do_nifti_image_metadata_elements_match("qoffset_x",       im1_sptr->qoffset_x,       im2_sptr->qoffset_x       ) &&
            do_nifti_image_metadata_elements_match("qoffset_y",       im1_sptr->qoffset_y,       im2_sptr->qoffset_y       ) &&
            do_nifti_image_metadata_elements_match("qoffset_z",       im1_sptr->qoffset_z,       im2_sptr->qoffset_z       ) &&
            do_nifti_image_metadata_elements_match("quatern_b",       im1_sptr->quatern_b,       im2_sptr->quatern_b       ) &&
            do_nifti_image_metadata_elements_match("quatern_c",       im1_sptr->quatern_c,       im2_sptr->quatern_c       ) &&
            do_nifti_image_metadata_elements_match("quatern_d",       im1_sptr->quatern_d,       im2_sptr->quatern_d       ) &&
            do_nifti_image_metadata_elements_match("scl_inter",       im1_sptr->scl_inter,       im2_sptr->scl_inter       ) &&
            do_nifti_image_metadata_elements_match("scl_slope",       im1_sptr->scl_slope,       im2_sptr->scl_slope       ) &&
            do_nifti_image_metadata_elements_match("sform_code",      im1_sptr->sform_code,      im2_sptr->sform_code      ) &&
            do_nifti_image_metadata_elements_match("slice_code",      im1_sptr->slice_code,      im2_sptr->slice_code      ) &&
            do_nifti_image_metadata_elements_match("slice_dim",       im1_sptr->slice_dim,       im2_sptr->slice_dim       ) &&
            do_nifti_image_metadata_elements_match("slice_duration",  im1_sptr->slice_duration,  im2_sptr->slice_duration  ) &&
            do_nifti_image_metadata_elements_match("slice_end",       im1_sptr->slice_end,       im2_sptr->slice_end       ) &&
            do_nifti_image_metadata_elements_match("slice_start",     im1_sptr->slice_start,     im2_sptr->slice_start     ) &&
            do_nifti_image_metadata_elements_match("swapsize",        im1_sptr->swapsize,        im2_sptr->swapsize        ) &&
            do_nifti_image_metadata_elements_match("time_units",      im1_sptr->time_units,      im2_sptr->time_units      ) &&
            do_nifti_image_metadata_elements_match("toffset",         im1_sptr->toffset,         im2_sptr->toffset         ) &&
            do_nifti_image_metadata_elements_match("xyz_units",       im1_sptr->xyz_units,       im2_sptr->xyz_units       ) &&
            do_nifti_image_metadata_elements_match("qto_ijk",         im1_sptr->qto_ijk,         im2_sptr->qto_ijk         ) &&
            do_nifti_image_metadata_elements_match("qto_xyz",         im1_sptr->qto_xyz,         im2_sptr->qto_xyz         ) &&
            do_nifti_image_metadata_elements_match("sto_ijk",         im1_sptr->sto_ijk,         im2_sptr->sto_ijk         ) &&
            do_nifti_image_metadata_elements_match("sto_xyz",         im1_sptr->sto_xyz,         im2_sptr->sto_xyz         );

    for (int i=0; i<8; i++) {
        if (!do_nifti_image_metadata_elements_match("dim["+std::to_string(i)+"]",    im1_sptr->dim[i],    im2_sptr->dim[i] ))   images_match = false;
        if (!do_nifti_image_metadata_elements_match("pixdim["+std::to_string(i)+"]", im1_sptr->pixdim[i], im2_sptr->pixdim[i])) images_match = false;
    }

#ifndef NDEBUG
    if (images_match) std::cout << "\tOK!\n";
#endif

    return images_match;
}

template<typename T>
bool NiftiImageData::do_nifti_image_metadata_elements_match(const std::string &name, const T &elem1, const T &elem2)
{
    if(float(fabs(elem1-elem2)) < 1.e-7F)
        return true;
    std::cout << "mismatch in " << name << " , (values: " <<  elem1 << " and " << elem2 << ")\n";
    return false;
}
template bool NiftiImageData::do_nifti_image_metadata_elements_match<float> (const std::string &name, const float &elem1, const float &elem2);

bool NiftiImageData::do_nifti_image_metadata_elements_match(const std::string &name, const mat44 &elem1, const mat44 &elem2)
{
    SIRFRegAffineTransformation e1(elem1.m), e2(elem2.m);
    if(e1 == e2)
        return true;
    std::cout << "mismatch in " << name << "\n";
    SIRFRegAffineTransformation::print({e1, e2});
    std::cout << "\n";
    return false;
}

/// Dump info of multiple nifti images
void NiftiImageData::dump_headers(const std::vector<NiftiImageData> &ims)
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
        images.push_back(ims[i].get_raw_nifti_sptr());

    // Print transformation matrices
    std::vector<SIRFRegAffineTransformation> qto_ijk_vec, qto_xyz_vec, sto_ijk_vec, sto_xyz_vec;
    for(unsigned j=0; j<images.size(); j++) {
        qto_ijk_vec.push_back(images[j]->qto_ijk.m);
        qto_xyz_vec.push_back(images[j]->qto_xyz.m);
        sto_ijk_vec.push_back(images[j]->sto_ijk.m);
        sto_xyz_vec.push_back(images[j]->sto_xyz.m);
    }
    std::cout << "\t" << std::left << std::setw(19) << "qto_ijk:" << "\n";
    SIRFRegAffineTransformation::print(qto_ijk_vec);
    std::cout << "\t" << std::left << std::setw(19) << "qto_xyz:" << "\n";
    SIRFRegAffineTransformation::print(qto_xyz_vec);
    std::cout << "\t" << std::left << std::setw(19) << "sto_ijk:" << "\n";
    SIRFRegAffineTransformation::print(sto_ijk_vec);
    std::cout << "\t" << std::left << std::setw(19) << "sto_xyz:" << "\n";
    SIRFRegAffineTransformation::print(sto_xyz_vec);

    // Print min
    std::string min_header = "min: ";
    std::cout << "\t" << std::left << std::setw(19) << min_header;
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i].get_min();

    // Print max
    std::cout << "\n\t" << std::left << std::setw(19) << "max: ";
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i].get_max();

    // Print mean
    std::cout << "\n\t" << std::left << std::setw(19) << "mean: ";
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i].get_mean();

    std::cout << "\n\n";
}

template<typename T>
void NiftiImageData::dump_nifti_element(const std::vector<NiftiImageData> &ims, const std::string &name, const T &call_back)
{
    std::string header = name + ": ";
    std::cout << "\t" << std::left << std::setw(19) << header;
    for(unsigned i=0; i<ims.size(); i++)
        std::cout << std::setw(19) << ims[i].get_raw_nifti_sptr().get()->*call_back;
    std::cout << "\n";
}

template<typename T>
void NiftiImageData::dump_nifti_element(const std::vector<NiftiImageData> &ims, const std::string &name, const T &call_back, const unsigned num_elems)
{
    for(unsigned i=0; i<num_elems; i++) {
        std::string header = name + "[" + std::to_string(i) + "]: ";
        std::cout << "\t" << std::left << std::setw(19) << header;
        for(unsigned j=0; j<ims.size(); j++)
            std::cout << std::setw(19) << (ims[j].get_raw_nifti_sptr().get()->*call_back)[i];
        std::cout << "\n";
    }
}

bool NiftiImageData::are_equal_to_given_accuracy(const NiftiImageData &im2, const float required_accuracy_compared_to_max) const
{
    const NiftiImageData &im1 = *this;

    if(!im1.is_initialised())
        throw std::runtime_error("NiftiImageData::are_equal_to_given_accuracy: Image 1 not initialised.");
    if(!im2.is_initialised())
        throw std::runtime_error("NiftiImageData::are_equal_to_given_accuracy: Image 2 not initialised.");

    // Get norm between two images
    float norm = im1.get_norm(im2);
    float epsilon = (im1.get_max()+im2.get_max())/2.F;
    epsilon *= required_accuracy_compared_to_max;

    if (norm < epsilon)
        return true;

    std::cout << "\nImages are not equal (norm > epsilon).\n";
    std::cout << "\tmax1                              = " << im1.get_max() << "\n";
    std::cout << "\tmax2                              = " << im1.get_max() << "\n";
    std::cout << "\tmin1                              = " << im1.get_min() << "\n";
    std::cout << "\tmin2                              = " << im2.get_min() << "\n";
    std::cout << "\trequired accuracy compared to max = " << required_accuracy_compared_to_max << "\n";
    std::cout << "\tepsilon                           = " << epsilon << "\n";
    std::cout << "\tnorm                              = " << norm << "\n";
    return false;
}
