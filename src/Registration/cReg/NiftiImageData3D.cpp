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

#include "sirf/cReg/NiftiImageData3D.h"
#include <nifti1_io.h>
#include <_reg_tools.h>
#include "sirf/cReg/SIRFRegAffineTransformation.h"

using namespace sirf;

template<class dataType>
NiftiImageData3D<dataType>::NiftiImageData3D(const ImageData& id)
{
    // Figure out dimensions
    Dimensions dim_map = id.dimensions();
    int dims[8];
    dims[0] = 3;
    dims[1] = dim_map["x"];
    dims[2] = dim_map["y"];
    dims[3] = dim_map["z"];
    dims[4] = 1;
    dims[5] = 1;
    dims[6] = 1;
    dims[7] = 1;

    // Create image
    nifti_image *im = nifti_make_new_nim(dims, DT_FLOAT32, 1);
    this->_nifti_image = std::shared_ptr<nifti_image>(im, nifti_image_free);

    // Spacing etc. are unknown until GeometricalInfo has been implemented

    // Copy data
    this->copy(id.begin(), this->begin(), this->end());

    // Check everything is ok
    reg_checkAndCorrectDimension(this->_nifti_image.get());

    // Always float
    this->set_up_data(NIFTI_TYPE_FLOAT32);
}

/* TODO UNCOMMENT WHEN GEOMETRICAL INFO IS IMPLEMENTED
template<class dataType>
bool NiftiImageData3D<dataType>::check_images_are_aligned(const VoxelisedGeometricalInfo3D &info) const
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
    if (this->_nifti_image->dim[1] != int(size[0])) ok_size = false;
    if (this->_nifti_image->dim[2] != int(size[1])) ok_size = false;
    if (this->_nifti_image->dim[3] != int(size[2])) ok_size = false;
    if (!ok_size) {
        std::cout << "\nWarning: Size does not match, can't fill image.\n";
        std::cout << "\tSTIR image   = (" << size[0]       << ", " << size[1]       << ", " << size[2]       << ")\n";
        std::cout << "\tNiftiImageData3D = (" << this->_nifti_image->dim[1] << ", " << this->_nifti_image->dim[2] << ", " << this->_nifti_image->dim[3] << ")\n";
    }

    // Check spacing
    bool ok_spacing = true;
    if (fabs(this->_nifti_image->dx - spacing[0]) > 1.e-7F) ok_spacing = false;
    if (fabs(this->_nifti_image->dy - spacing[1]) > 1.e-7F) ok_spacing = false;
    if (fabs(this->_nifti_image->dz - spacing[2]) > 1.e-7F) ok_spacing = false;
    if (!ok_spacing) {
        std::cout << "\nWarning: Spacing does not match, can't fill image.\n";
        std::cout << "\tSTIR image   = (" << spacing[0]       << ", " << spacing[1]       << ", " << spacing[2]       << ")\n";
        std::cout << "\tNiftiImageData3D = (" << this->_nifti_image->dx << ", " << this->_nifti_image->dy << ", " << this->_nifti_image->dz << ")\n";
    }

    // Check qto_xyz
    SIRFRegAffineTransformation<dataType> qto_xyz(this->_nifti_image->qto_xyz.m);
    SIRFRegAffineTransformation<dataType> stir_qto_xyz;
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            stir_qto_xyz[i][j] = tm[i][j];
    bool ok_qto_xyz = (stir_qto_xyz == qto_xyz);
    if (!ok_qto_xyz) {
        std::cout << "\nWarning: qto_xyz does not match, can't fill image.\n";
        std::vector<SIRFRegAffineTransformation<dataType> > mats = {stir_qto_xyz, qto_xyz};
        SIRFRegAffineTransformation<dataType>::print(mats);
    }

    // Return if everything is ok or not.
    return (ok_size && ok_spacing && ok_qto_xyz);
}*/

namespace sirf {
template class NiftiImageData3D<float>;
}
