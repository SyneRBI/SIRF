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

#include "NiftiImage3D.h"
#include "SIRFRegMat44.h"
#include "SIRFRegMisc.h"
#include <nifti1_io.h>
#include <_reg_tools.h>
#include "stir_data_containers.h"

using namespace sirf;

NiftiImage3D::NiftiImage3D(const PETImageData &pet_image)
{
    std::cout << "Converting PET image to nifti image..." << std::flush;

    // Set up the nifti
    set_up_nifti(pet_image.get_patient_coord_geometrical_info());

    // Copy the data (both datatypes are float)
    pet_image.get_data(_data);

    std::cout << "Done!\n";
}

NiftiImage3D::NiftiImage3D(const MRImageData &)
{
    std::cout << "\n\nTODO\n\n";
    exit(0);
}

void NiftiImage3D::set_up_nifti(const VoxelisedGeometricalInfo3D &info)
{
    typedef VoxelisedGeometricalInfo3D Info;
    Info::Size            size    = info.get_size();
    Info::Spacing         spacing = info.get_spacing();
    Info::TransformMatrix tm      = info.calculate_index_to_physical_point_matrix();

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
    _nifti_image = std::shared_ptr<nifti_image>(im, nifti_image_free);

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
    _nifti_image->pixdim[0]=_nifti_image->qfac;

    // Check everything is ok
    reg_checkAndCorrectDimension(_nifti_image.get());

    // Always float
    set_up_data(NIFTI_TYPE_FLOAT32);
}

void NiftiImage3D::copy_data_to(PETImageData &pet_image) const
{
    std::cout << "Filling PET image from nifti image..." << std::flush;

    if (!check_images_are_aligned(
                pet_image.get_patient_coord_geometrical_info()))
        throw std::runtime_error("Cannot copy data from NiftImage to STIRImageData as they are not aligned.");

    pet_image.set_data(_data);

    std::cout << "Done!\n";
}

void NiftiImage3D::copy_data_to(MRImageData &) const
{
    std::cout << "\n\nTODO\n\n";
    exit(0);
}

bool NiftiImage3D::check_images_are_aligned(const VoxelisedGeometricalInfo3D &info) const
{
    // Check the nifti exists
    if (!this->is_initialised()) {
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
    if (_nifti_image->dim[1] != int(size[0])) ok_size = false;
    if (_nifti_image->dim[2] != int(size[1])) ok_size = false;
    if (_nifti_image->dim[3] != int(size[2])) ok_size = false;
    if (!ok_size) {
        std::cout << "\nWarning: Size does not match, can't fill image.\n";
        std::cout << "\tSTIR image   = (" << size[0]       << ", " << size[1]       << ", " << size[2]       << ")\n";
        std::cout << "\tNiftiImage3D = (" << _nifti_image->dim[1] << ", " << _nifti_image->dim[2] << ", " << _nifti_image->dim[3] << ")\n";
    }

    // Check spacing
    bool ok_spacing = true;
    if (fabs(_nifti_image->dx - spacing[0]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->dy - spacing[1]) > 1.e-7F) ok_spacing = false;
    if (fabs(_nifti_image->dz - spacing[2]) > 1.e-7F) ok_spacing = false;
    if (!ok_spacing) {
        std::cout << "\nWarning: Spacing does not match, can't fill image.\n";
        std::cout << "\tSTIR image   = (" << spacing[0]       << ", " << spacing[1]       << ", " << spacing[2]       << ")\n";
        std::cout << "\tNiftiImage3D = (" << _nifti_image->dx << ", " << _nifti_image->dy << ", " << _nifti_image->dz << ")\n";
    }

    // Check qto_xyz
    SIRFRegMat44 qto_xyz(_nifti_image->qto_xyz);
    SIRFRegMat44 stir_qto_xyz;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            stir_qto_xyz[i][j] = tm[i][j];
    bool ok_qto_xyz = (stir_qto_xyz == qto_xyz);
    if (!ok_qto_xyz) {
        std::cout << "\nWarning: qto_xyz does not match, can't fill image.\n";
        std::vector<SIRFRegMat44> mats = {stir_qto_xyz, qto_xyz};
        SIRFRegMat44::print(mats);
    }

    // Return if everything is ok or not.
    return (ok_size && ok_spacing && ok_qto_xyz);
}
