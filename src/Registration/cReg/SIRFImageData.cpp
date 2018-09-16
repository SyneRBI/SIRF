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

#include "SIRFImageData.h"
#include "SIRFRegMisc.h"
#include <nifti1_io.h>
#include <_reg_tools.h>
#include "stir_data_containers.h"

using namespace std;
using namespace sirf;

SIRFImageData SIRFImageData::operator=(const SIRFImageData& to_copy)
{
    // Check for self-assignment
    if (this != &to_copy)
        SIRFRegMisc::copy_nifti_image(_nifti_image,to_copy._nifti_image);

    return *this;
}

SIRFImageData::SIRFImageData(const std::string &filename)
{
    SIRFRegMisc::open_nifti_image(_nifti_image,filename);
}

SIRFImageData::SIRFImageData(const nifti_image &image_nifti)
{
    SIRFRegMisc::copy_nifti_image(_nifti_image,make_shared<nifti_image>(image_nifti));
    reg_checkAndCorrectDimension(_nifti_image.get());
}

SIRFImageData::SIRFImageData(const std::shared_ptr<nifti_image> image_nifti)
{
    SIRFRegMisc::copy_nifti_image(_nifti_image,image_nifti);
    reg_checkAndCorrectDimension(_nifti_image.get());
}

SIRFImageData::SIRFImageData(const sirf::PETImageData &pet_image)
{
    cout << "Converting PET image to nifti image..." << flush;

    // Set up the nifti
    set_up_nifti(pet_image.get_patient_coord_geometrical_info());

    // Copy the data. this cast is ok because PETImageData is always float.
    float *data = static_cast<float*>(_nifti_image->data);
    pet_image.get_data(data);

    cout << "Done!\n";
}

SIRFImageData::SIRFImageData(const MRImageData &)
{
    cout << "\n\nTODO\n\n";
    exit(0);
}

SIRFImageData SIRFImageData::operator+ (const SIRFImageData& c) const
{
    if (!this->is_initialised())
        throw runtime_error("Can't add SIRFImageData as first image is not initialised.");
    if (!c.is_initialised())
        throw runtime_error("Can't add SIRFImageData as second image is not initialised.");

    if (!SIRFRegMisc::do_nifti_image_metadata_match(*this, c))
        throw runtime_error("Can't add SIRFImageData as metadata do not match.");

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
    ss << "SIRFImageData::operator+ not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

SIRFImageData SIRFImageData::operator- (const SIRFImageData& c) const
{
    if (!this->is_initialised())
        throw runtime_error("Can't subtract SIRFImageData as first image is not initialised.");
    if (!c.is_initialised())
        throw runtime_error("Can't subtract SIRFImageData as second image is not initialised.");

    if (!SIRFRegMisc::do_nifti_image_metadata_match(*this, c))
        throw runtime_error("Can't subtract SIRFImageData as metadata do not match.");

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
    ss << "SIRFImageData::operator- not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

void SIRFImageData::set_up_nifti(const VoxelisedGeometricalInfo3D &info)
{
    typedef VoxelisedGeometricalInfo3D Info;
    Info::Size            size    = info.get_size();
    Info::Spacing         spacing = info.get_spacing();
    Info::TransformMatrix tm      = info.calculate_index_to_physical_point_matrix();

    _nifti_image = std::shared_ptr<nifti_image>(new nifti_image, nifti_image_free);
    _nifti_image->dim[0]=_nifti_image->ndim=3;
    // Size
    _nifti_image->dim[1]=_nifti_image->nx=int(size[0]);
    _nifti_image->dim[2]=_nifti_image->ny=int(size[1]);
    _nifti_image->dim[3]=_nifti_image->nz=int(size[2]);
    _nifti_image->dim[4]=_nifti_image->nt=1;
    _nifti_image->dim[5]=_nifti_image->nu=1;
    _nifti_image->dim[6]=_nifti_image->nv=1;
    _nifti_image->dim[7]=_nifti_image->nw=1;
    _nifti_image->nvox=unsigned(_nifti_image->nx*_nifti_image->ny*_nifti_image->nz*_nifti_image->nt*_nifti_image->nu);
    // Spacing (extra dimensions are 0 by default)
    _nifti_image->pixdim[1]=_nifti_image->dx=spacing[0];
    _nifti_image->pixdim[2]=_nifti_image->dy=spacing[1];
    _nifti_image->pixdim[3]=_nifti_image->dz=spacing[2];
    // Data types
    _nifti_image->datatype = DT_FLOAT32;
    _nifti_image->nbyper = sizeof(float);
    _nifti_image->swapsize = sizeof(float);
    _nifti_image->intent_code = NIFTI_INTENT_NONE;
    _nifti_image->xyz_units=2; // distances in mm
    _nifti_image->nifti_type=1;
    _nifti_image->byteorder=1;
    _nifti_image->scl_inter=0.F;
    _nifti_image->scl_slope=1.F;
    _nifti_image->iname_offset=352;
    // Set the transformation matrix information
    _nifti_image->qform_code=1;
    for (int i=0;i<4;++i)
        for (int j=0;j<4;++j)
            _nifti_image->qto_xyz.m[i][j]=tm[i][j];
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
    // Null pointers for some stuff
    _nifti_image->fname = nullptr;
    _nifti_image->iname = nullptr;
    _nifti_image->num_ext = 0;
    _nifti_image->ext_list = nullptr;

    // Check everything is ok
    reg_checkAndCorrectDimension(_nifti_image.get());

    // Allocate the data
    _nifti_image->data = static_cast<void *>(calloc(_nifti_image->nvox, unsigned(_nifti_image->nbyper)));
}

std::shared_ptr<nifti_image> SIRFImageData::get_raw_nifti_sptr() const
{
    if (!this->is_initialised())
        throw runtime_error("Warning, nifti has not been initialised.");
    return _nifti_image;
}

void SIRFImageData::copy_data_to(sirf::PETImageData &pet_image) const
{
    cout << "Filling PET image from nifti image..." << flush;

    bool everything_ok =
            check_images_are_aligned(
                pet_image.get_patient_coord_geometrical_info());

    if (!everything_ok)
        throw std::runtime_error("Cannot copy data from SIRFImageData to STIRImageData as they are not aligned.");

    // If datatype is already float
    if (_nifti_image->datatype == DT_FLOAT32) {
        float *nifti_data_ptr = static_cast<float *>(_nifti_image->data);
        pet_image.set_data(nifti_data_ptr);
    }
    // If not, cast to it
    else {
        SIRFImageData temp = this->deep_copy();
        SIRFRegMisc::change_datatype<float>(temp);
        float *nifti_data_ptr = static_cast<float *>(temp.get_raw_nifti_sptr()->data);
        pet_image.set_data(nifti_data_ptr);
    }

    cout << "Done!\n";
}

void SIRFImageData::copy_data_to(MRImageData &) const
{
    cout << "\n\nTODO\n\n";
    exit(0);
}

bool SIRFImageData::check_images_are_aligned(const VoxelisedGeometricalInfo3D &info) const
{
    // Check the nifti exists
    if (!_nifti_image) {
        std::cout << "\nWarning: Nifti image not initialised, can't fill image.\n";
        return false;
    }

    // Check the info all matches (they should have resampled first)
    typedef VoxelisedGeometricalInfo3D Info;
    Info::Size            size    = info.get_size();
    Info::Spacing         spacing = info.get_spacing();
    Info::TransformMatrix tm      = info.calculate_index_to_physical_point_matrix();

    // Check size
    bool ok_size = true;
    if (_nifti_image->dim[0] != 3)                       ok_size = false;
    if (_nifti_image->dim[1] != int(size[0]))            ok_size = false;
    if (_nifti_image->dim[2] != int(size[1]))            ok_size = false;
    if (_nifti_image->dim[3] != int(size[2]))            ok_size = false;
    if (_nifti_image->dim[4] != 1)                       ok_size = false;
    if (_nifti_image->dim[5] != 1)                       ok_size = false;
    if (_nifti_image->dim[6] != 1)                       ok_size = false;
    if (_nifti_image->dim[7] != 1)                       ok_size = false;
    if (_nifti_image->nx     != int(size[0]))            ok_size = false;
    if (_nifti_image->ny     != int(size[1]))            ok_size = false;
    if (_nifti_image->nz     != int(size[2]))            ok_size = false;
    if (_nifti_image->nt     != 1)                       ok_size = false;
    if (_nifti_image->nu     != 1)                       ok_size = false;
    if (_nifti_image->nv     != 1)                       ok_size = false;
    if (_nifti_image->nw     != 1)                       ok_size = false;
    if (_nifti_image->nvox   != size[0]*size[1]*size[2]) ok_size = false;
    if (!ok_size)
        std::cout << "\nWarning: Size does not match, can't fill image.\n";

    // Check spacing
    bool ok_spacing = true;
    if (fabs(_nifti_image->pixdim[1] - spacing[0]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->pixdim[2] - spacing[1]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->pixdim[3] - spacing[2]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->dx        - spacing[0]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->dy        - spacing[1]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->dz        - spacing[2]) > 1.e-7F) ok_spacing = false;
    if (!ok_spacing)
        std::cout << "\nWarning: Spacing does not match, can't fill image.\n";

    // Check offsets
    bool ok_offset = true;
    if (fabs(tm[0][3] - _nifti_image->qoffset_x) > 1.e-7F) ok_offset = false;
    if (fabs(tm[1][3] - _nifti_image->qoffset_y) > 1.e-7F) ok_offset = false;
    if (fabs(tm[2][3] - _nifti_image->qoffset_z) > 1.e-7F) ok_offset = false;
    if (!ok_offset)
        std::cout << "\nWarning: qoffset does not match, can't fill image.\n";

    // Check qto_xyz
    bool ok_qto_xyz = true;
    for (int i=0;i<4;++i)
        for (int j=0;j<4;++j)
            if (fabs(_nifti_image->qto_xyz.m[i][j] - tm[i][j]) > 1.e-7F)
                ok_qto_xyz = false;
    if (!ok_qto_xyz)
        std::cout << "\nWarning: qto_xyz does not match, can't fill image.\n";

    // Check qto_ijk
    bool ok_qto_ijk = true;
    if (fabs( _nifti_image->qto_ijk.m[0][0] - 1.F/tm[0][0])        > 1.e-7F) ok_qto_ijk = false;
    if (fabs( _nifti_image->qto_ijk.m[1][1] - 1.F/tm[1][1])        > 1.e-7F) ok_qto_ijk = false;
    if (fabs( _nifti_image->qto_ijk.m[2][2] - 1.F/tm[2][2])        > 1.e-7F) ok_qto_ijk = false;
    if (fabs( _nifti_image->qto_ijk.m[0][3] - tm[0][3]/spacing[0]) > 1.e-7F) ok_qto_ijk = false;
    if (fabs( _nifti_image->qto_ijk.m[1][3] - tm[1][3]/spacing[1]) > 1.e-7F) ok_qto_ijk = false;
    if (fabs( _nifti_image->qto_ijk.m[2][3] - tm[2][3]/spacing[2]) > 1.e-7F) ok_qto_ijk = false;
    if (!ok_qto_ijk)
        std::cout << "\nWarning: qto_ijk does not match, can't fill image.\n";

    // Check datatype float
    bool ok_datatype = true;
    if (_nifti_image->datatype != DT_FLOAT32   ) ok_datatype = false;
    if ( _nifti_image->nbyper  != sizeof(float)) ok_datatype = false;
    if (!ok_datatype)
        std::cout << "\nWarning: datatype is not float, can't fill image.\n";

    // Return if not everything is ok
    if(ok_size && ok_spacing && ok_offset && ok_qto_xyz && ok_qto_ijk && ok_datatype)
        return true;
    return false;
}

void SIRFImageData::save_to_file(const std::string &filename) const
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

float SIRFImageData::get_max() const
{
    if(!_nifti_image)
        throw runtime_error("Image not initialised.");

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
    ss << "SIRFImageData::get_max not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

float SIRFImageData::get_min() const
{
    if(!_nifti_image)
        throw runtime_error("Image not initialised.");

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
    ss << "SIRFImageData::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

float SIRFImageData::get_element(int x, int y, int z, int t, int u, int v, int w) const
{
    if(!this->is_initialised())
        throw runtime_error("Image not initialised.");

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
    ss << "SIRFImageData::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

float SIRFImageData::get_sum() const
{
    if(!_nifti_image)
        throw runtime_error("Image not initialised.");

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
    ss << "SIRFImageData::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

void SIRFImageData::fill(const float &v)
{
    if(!this->is_initialised())
        throw runtime_error("Image not initialised.");

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
    ss << "SIRFImageData::get_min not implemented for your data type: ";
    ss << nifti_datatype_string(_nifti_image->datatype);
    ss << " (bytes per voxel: ";
    ss << _nifti_image->nbyper << ").";
    throw std::runtime_error(ss.str());
}

SIRFImageData SIRFImageData::deep_copy() const
{
    SIRFImageData copy;
    copy = *this;
    return copy;
}
