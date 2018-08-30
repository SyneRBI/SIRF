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

#include "SIRFImageDataDeformation.h"
#include "SIRFRegMisc.h"

using namespace std;

SIRFImageDataDeformation::SIRFImageDataDeformation(const std::string &filename)
    : SIRFImageData(filename)
{
    if (_nifti_image->ndim != 5)
        throw runtime_error("A deformation/displacement field image should have ndim=5.");
}

SIRFImageDataDeformation::SIRFImageDataDeformation(const nifti_image *image_nifti)
    : SIRFImageData(image_nifti)
{
    if (_nifti_image->ndim != 5)
        throw runtime_error("A deformation/displacement field image should have ndim=5.");
}

SIRFImageDataDeformation::SIRFImageDataDeformation(const std::shared_ptr<nifti_image> image_nifti)
    : SIRFImageData(image_nifti)
{
    if (_nifti_image->ndim != 5)
        throw runtime_error("A deformation/displacement field image should have ndim=5.");
}

void SIRFImageDataDeformation::create_from_3D_image(const SIRFImageData &image)
{
    std::shared_ptr<nifti_image> image_sptr = image.get_image_as_nifti();

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

    _nifti_image = make_shared<nifti_image>(*output_ptr);
}

void SIRFImageDataDeformation::save_to_file(const std::string &filename, bool split_xyz, string type)
{
    // Check that the disp image exists
    if (!_nifti_image)
        throw std::runtime_error("Error, " + type + " image not available. Have you run the registration?");

    // Check that filename isn't blank
    if (filename == "")
        throw std::runtime_error("Error, cannot write " + type + " image to file because filename is blank.");

    cout << "\nSaving " + type + " image to file (" << filename << ")..." << flush;

    SIRFRegMisc::save_multicomponent_nifti_image(_nifti_image,filename,split_xyz);

    cout << "Done.\n";
}
