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
\brief Class for deformation and displacement SIRF image data.

\author Richard Brown
\author CCP PETMR
*/

#include "NiftiImage3DTensor.h"
#include "SIRFRegMisc.h"
#include <boost/format.hpp>

using namespace std;
using namespace sirf;

NiftiImage3DTensor::NiftiImage3DTensor(const NiftiImage3D &x, const NiftiImage3D &y, const NiftiImage3D &z)
{
    // Check everything is intialised
    if (!x.is_initialised() || !y.is_initialised() || !z.is_initialised())
        throw runtime_error("NiftiImage3DTensor: x,y,z->tensor: Can't create from separate 3D components, as some are uninitialised.");

    if (!SIRFRegMisc::do_nifti_image_metadata_match(x,y))
        throw runtime_error("NiftiImage3DTensor: x,y,z->tensor: x and y components don't match.");
    if (!SIRFRegMisc::do_nifti_image_metadata_match(x,z))
        throw runtime_error("NiftiImage3DTensor: x,y,z->tensor: x and z components don't match.");

    // Create a 4D from one of the components
    this->create_from_3D_image(x);
    vector<NiftiImage3D> ims{x, y, z};

    // for nu=3, the tensor data is stored last.
    //So memcpy x into first third, y into second third and z into last third
    size_t mem = x.get_raw_nifti_sptr()->nvox*size_t(x.get_raw_nifti_sptr()->nbyper);
    char *dest = static_cast<char*>(_nifti_image->data);
    for (size_t i=0; i<3; ++i) {
        size_t index = mem*i;
        char *src  = static_cast<char*>(ims[i].get_raw_nifti_sptr()->data);
        memcpy(dest+index, src, mem);
    }
}

void NiftiImage3DTensor::create_from_3D_image(const NiftiImage3D &image)
{
    if (!image.is_initialised())
        throw runtime_error("NiftiImage3DTensor::create_from_3D_image. Input image not initialised.");

    std::shared_ptr<nifti_image> image_sptr = image.get_raw_nifti_sptr();

    // Calculate deformation field image
    nifti_image *output_ptr;
    output_ptr = nifti_copy_nim_info(image_sptr.get());
    output_ptr->dim[0]=output_ptr->ndim=5;
    output_ptr->dim[1]=output_ptr->nx=image_sptr->nx;
    output_ptr->dim[2]=output_ptr->ny=image_sptr->ny;
    output_ptr->dim[3]=output_ptr->nz=image_sptr->nz;
    output_ptr->dim[4]=output_ptr->nt=1;output_ptr->pixdim[4]=output_ptr->dt=1.0;
    if(image_sptr->nz>1) output_ptr->dim[5]=output_ptr->nu=3;
    else output_ptr->dim[5]=output_ptr->nu=2;
    output_ptr->pixdim[5]=output_ptr->du=1.0;
    output_ptr->dim[6]=output_ptr->nv=1;output_ptr->pixdim[6]=output_ptr->dv=1.0;
    output_ptr->dim[7]=output_ptr->nw=1;output_ptr->pixdim[7]=output_ptr->dw=1.0;
    output_ptr->nvox=output_ptr->nx*output_ptr->ny*output_ptr->nz*output_ptr->nt*output_ptr->nu;
    output_ptr->datatype = DT_FLOAT32;
    output_ptr->nbyper = sizeof(float);
    output_ptr->data = static_cast<void *>(calloc(output_ptr->nvox, unsigned(output_ptr->nbyper)));
    output_ptr->intent_code = NIFTI_INTENT_VECTOR;

    _nifti_image = std::shared_ptr<nifti_image>(output_ptr, nifti_image_free);
}

void NiftiImage3DTensor::save_to_file_split_xyz_components(const std::string &filename) const
{
    // Check that the disp image exists
    if (!this->is_initialised())
        throw std::runtime_error("Error, cannot write " + filename + " to file because it has not been initialised.");

    // Check that filename isn't blank
    if (filename.size() == 0)
        throw std::runtime_error("Error, cannot write " + filename + " to file because filename is blank.");

    cout << "\nSaving to file (" << filename << ")..." << flush;

    SIRFRegMisc::save_multicomponent_nifti_image_split_xyz(_nifti_image,filename);

    cout << "Done.\n";
}

void NiftiImage3DTensor::save_to_file_split_xyz_components(const std::string &filename_x, const std::string &filename_y, const std::string &filename_z) const
{
    // Check that the disp image exists
    if (!this->is_initialised())
        throw std::runtime_error("Error, cannot write to " + filename_x + ", " + filename_y + ", " + filename_z + " because it has not been initialised.");

    // Check that filename isn't blank
    if (filename_x.size() == 0 || filename_y.size() == 0 || filename_z.size() == 0)
        throw std::runtime_error("Error, cannot write to " + filename_x + ", " + filename_y + ", " + filename_z + " because one or more filenames is blank.");

    cout << "\nSaving to file (" << filename_x + ", " + filename_y + ", " + filename_z << ")..." << flush;

    SIRFRegMisc::save_multicomponent_nifti_image_split_xyz(_nifti_image,filename_x,filename_y,filename_z);

    cout << "Done.\n";
}
