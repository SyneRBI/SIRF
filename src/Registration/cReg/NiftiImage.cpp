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

#include "NiftiImage.h"
#include "SIRFRegMisc.h"
#include <nifti1_io.h>
#include <_reg_tools.h>
#include "stir_data_containers.h"
#include "NiftiImage3DTensor.h"
#include "NiftiImage3DDeformation.h"
#include "NiftiImage3DDisplacement.h"

using namespace std;
using namespace sirf;

NiftiImage NiftiImage::operator=(const NiftiImage& to_copy)
{
    // Check for self-assignment
    if (this != &to_copy) {
        if (to_copy.is_initialised())
            SIRFRegMisc::copy_nifti_image(_nifti_image,to_copy._nifti_image);
        else
            throw runtime_error("trying to copy empty image");
    }

    return *this;
}

NiftiImage::NiftiImage(const std::string &filename)
{
    SIRFRegMisc::open_nifti_image(_nifti_image,filename);
}

NiftiImage::NiftiImage(const nifti_image &image_nifti)
{
    SIRFRegMisc::copy_nifti_image(_nifti_image,make_shared<nifti_image>(image_nifti));
    reg_checkAndCorrectDimension(_nifti_image.get());
}

NiftiImage::NiftiImage(const std::shared_ptr<nifti_image> image_nifti)
{
    SIRFRegMisc::copy_nifti_image(_nifti_image,image_nifti);
    reg_checkAndCorrectDimension(_nifti_image.get());
}

NiftiImage NiftiImage::operator+ (const NiftiImage& c) const
{
    if (!this->is_initialised())
        throw runtime_error("Can't add NiftImage as first image is not initialised.");
    if (!c.is_initialised())
        throw runtime_error("Can't add NiftImage as second image is not initialised.");

    if (!SIRFRegMisc::do_nifti_image_metadata_match(*this, c))
        throw runtime_error("Can't add NiftImage as metadata do not match.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::sum_arrays<bool>              (*this,c);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::sum_arrays<signed char>       (*this,c);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::sum_arrays<signed short>      (*this,c);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::sum_arrays<signed int>        (*this,c);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::sum_arrays<float>             (*this,c);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::sum_arrays<double>            (*this,c);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::sum_arrays<unsigned char>     (*this,c);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::sum_arrays<unsigned short>    (*this,c);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::sum_arrays<unsigned int>      (*this,c);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::sum_arrays<signed long long>  (*this,c);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::sum_arrays<unsigned long long>(*this,c);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::sum_arrays<long double>       (*this,c);

    stringstream ss;
    ss << "NiftImage::operator+ not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

NiftiImage NiftiImage::operator- (const NiftiImage& c) const
{
    if (!this->is_initialised())
        throw runtime_error("Can't subtract NiftImage as first image is not initialised.");
    if (!c.is_initialised())
        throw runtime_error("Can't subtract NiftImage as second image is not initialised.");

    if (!SIRFRegMisc::do_nifti_image_metadata_match(*this, c))
        throw runtime_error("Can't subtract NiftImage as metadata do not match.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::sub_arrays<bool>              (*this,c);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::sub_arrays<signed char>       (*this,c);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::sub_arrays<signed short>      (*this,c);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::sub_arrays<signed int>        (*this,c);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::sub_arrays<float>             (*this,c);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::sub_arrays<double>            (*this,c);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::sub_arrays<unsigned char>     (*this,c);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::sub_arrays<unsigned short>    (*this,c);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::sub_arrays<unsigned int>      (*this,c);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::sub_arrays<signed long long>  (*this,c);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::sub_arrays<unsigned long long>(*this,c);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::sub_arrays<long double>       (*this,c);

    stringstream ss;
    ss << "NiftImage::operator- not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

std::shared_ptr<nifti_image> NiftiImage::get_raw_nifti_sptr() const
{
    if (!this->is_initialised())
        throw runtime_error("Warning, nifti has not been initialised.");
    return _nifti_image;
}

void NiftiImage::save_to_file(const std::string &filename) const
{
    if (!this->is_initialised())
        throw runtime_error("Cannot save image to file.");

    cout << "\nSaving image to file (" << filename << ")..." << flush;

    boost::filesystem::path filename_boost(filename);

    // If the folder doesn't exist, create it
    if (!boost::filesystem::exists(filename_boost.parent_path())) {
        if (filename_boost.parent_path().string() != "") {
            cout << "\n\tCreating folder: \"" << filename_boost.parent_path().string() << "\"\n" << flush;
            boost::filesystem::create_directory(filename_boost.parent_path());
        }
    }

    nifti_set_filenames(_nifti_image.get(), filename.c_str(), 0, 0);
    nifti_image_write(_nifti_image.get());
    cout << "done.\n\n";
}

float NiftiImage::get_max() const
{
    if(!this->is_initialised())
        throw runtime_error("NiftiImage::get_max(): Image not initialised.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::get_array_max<bool>              (*this);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::get_array_max<signed char>       (*this);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::get_array_max<signed short>      (*this);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::get_array_max<signed int>        (*this);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::get_array_max<float>             (*this);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::get_array_max<double>            (*this);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::get_array_max<unsigned char>     (*this);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::get_array_max<unsigned short>    (*this);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::get_array_max<unsigned int>      (*this);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::get_array_max<signed long long>  (*this);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::get_array_max<unsigned long long>(*this);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::get_array_max<long double>       (*this);

    stringstream ss;
    ss << "NiftImage::get_max not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

float NiftiImage::get_min() const
{
    if(!this->is_initialised())
        throw runtime_error("NiftiImage::get_min(): Image not initialised.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::get_array_min<bool>              (*this);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::get_array_min<signed char>       (*this);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::get_array_min<signed short>      (*this);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::get_array_min<signed int>        (*this);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::get_array_min<float>             (*this);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::get_array_min<double>            (*this);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::get_array_min<unsigned char>     (*this);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::get_array_min<unsigned short>    (*this);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::get_array_min<unsigned int>      (*this);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::get_array_min<signed long long>  (*this);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::get_array_min<unsigned long long>(*this);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::get_array_min<long double>       (*this);

    stringstream ss;
    ss << "NiftImage::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

float NiftiImage::get_element(int x, int y, int z, int t, int u, int v, int w) const
{
    if(!this->is_initialised())
        throw runtime_error("NiftiImage::get_element(): Image not initialised.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::get_array_element<bool>              (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::get_array_element<signed char>       (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::get_array_element<signed short>      (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::get_array_element<signed int>        (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::get_array_element<float>             (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::get_array_element<double>            (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::get_array_element<unsigned char>     (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::get_array_element<unsigned short>    (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::get_array_element<unsigned int>      (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::get_array_element<signed long long>  (*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::get_array_element<unsigned long long>(*this, x, y, z, t, u, v, w);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::get_array_element<long double>       (*this, x, y, z, t, u, v, w);

    stringstream ss;
    ss << "NiftImage::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

float NiftiImage::get_sum() const
{
    if(!this->is_initialised())
        throw runtime_error("NiftiImage::get_sum(): Image not initialised.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::get_array_sum<bool>              (*this);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::get_array_sum<signed char>       (*this);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::get_array_sum<signed short>      (*this);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::get_array_sum<signed int>        (*this);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::get_array_sum<float>             (*this);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::get_array_sum<double>            (*this);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::get_array_sum<unsigned char>     (*this);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::get_array_sum<unsigned short>    (*this);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::get_array_sum<unsigned int>      (*this);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::get_array_sum<signed long long>  (*this);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::get_array_sum<unsigned long long>(*this);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::get_array_sum<long double>       (*this);

    stringstream ss;
    ss << "NiftImage::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

void NiftiImage::fill(const float &v)
{
    if(!this->is_initialised())
        throw runtime_error("NiftiImage::fill(): Image not initialised.");

    if (_nifti_image->datatype == DT_BINARY)   return SIRFRegMisc::fill_array<bool>              (*this, v);
    if (_nifti_image->datatype == DT_INT8)     return SIRFRegMisc::fill_array<signed char>       (*this, v);
    if (_nifti_image->datatype == DT_INT16)    return SIRFRegMisc::fill_array<signed short>      (*this, v);
    if (_nifti_image->datatype == DT_INT32)    return SIRFRegMisc::fill_array<signed int>        (*this, v);
    if (_nifti_image->datatype == DT_FLOAT32)  return SIRFRegMisc::fill_array<float>             (*this, v);
    if (_nifti_image->datatype == DT_FLOAT64)  return SIRFRegMisc::fill_array<double>            (*this, v);
    if (_nifti_image->datatype == DT_UINT8)    return SIRFRegMisc::fill_array<unsigned char>     (*this, v);
    if (_nifti_image->datatype == DT_UINT16)   return SIRFRegMisc::fill_array<unsigned short>    (*this, v);
    if (_nifti_image->datatype == DT_UINT32)   return SIRFRegMisc::fill_array<unsigned int>      (*this, v);
    if (_nifti_image->datatype == DT_INT64)    return SIRFRegMisc::fill_array<signed long long>  (*this, v);
    if (_nifti_image->datatype == DT_UINT64)   return SIRFRegMisc::fill_array<unsigned long long>(*this, v);
    if (_nifti_image->datatype == DT_FLOAT128) return SIRFRegMisc::fill_array<long double>       (*this, v);

    stringstream ss;
    ss << "NiftImage::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

NiftiImage NiftiImage::deep_copy() const
{
    NiftiImage copy;
    copy = *this;
    return copy;
}

void NiftiImage::get_dimensions(int dims[8]) const
{
    for (int i=0; i<8; ++i)
        dims[i] = _nifti_image->dim[i];
}

void NiftiImage::check_dimensions(const NiftiImageType image_type)
{
    if (!this->is_initialised())
        throw runtime_error("NiftiImage::check_dimensions(): Image not initialised.");

    int ndim, nt, nu, intent_code, intent_p1;
    if        (image_type == _general)  { ndim=-1; nt=-1; nu=-1; intent_code = NIFTI_INTENT_NONE;   intent_p1=-1;         }
    else   if (image_type == _3D)       { ndim= 3; nt= 1; nu= 1; intent_code = NIFTI_INTENT_NONE;   intent_p1=-1;         }
    else   if (image_type == _3DTensor) { ndim= 5; nt= 1; nu= 3; intent_code = NIFTI_INTENT_VECTOR; intent_p1=-1;         }
    else   if (image_type == _3DDisp)   { ndim= 5; nt= 1; nu= 3; intent_code = NIFTI_INTENT_VECTOR; intent_p1=DISP_FIELD; }
    else /*if (image_type == _3DDef)*/  { ndim= 5; nt= 1; nu= 3; intent_code = NIFTI_INTENT_VECTOR; intent_p1=DEF_FIELD;  }

    // Check everthing is as it should be. -1 means we don't care about it
    // (e.g., NiftiImage3D doesn't care about intent_p1, which is used by NiftyReg for Disp/Def fields)
    bool everything_ok = true;
    if (ndim         != -1 && ndim        != _nifti_image->ndim)        everything_ok = false;
    if ( nu          != -1 && nu          != _nifti_image->nu  )        everything_ok = false;
    if ( nt          != -1 && nt          != _nifti_image->nt  )        everything_ok = false;
    if ( intent_code != -1 && intent_code != _nifti_image->intent_code) everything_ok = false;
    if ( intent_p1   != -1 && intent_p1   != _nifti_image->intent_p1)   everything_ok = false;

    if (everything_ok)
        return;

    // If not, throw an error.
    stringstream ss;
    ss << "Trying to construct a ";
    if      (typeid(*this) == typeid(NiftiImage3D))             ss << "NiftiImage3D";
    else if (typeid(*this) == typeid(NiftiImage3DTensor))       ss << "NiftiImage3DTensor";
    else if (typeid(*this) == typeid(NiftiImage3DDisplacement)) ss << "NiftiImage3DDisplacement";
    else if (typeid(*this) == typeid(NiftiImage3DDeformation))  ss << "NiftiImage3DDeformation";
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
