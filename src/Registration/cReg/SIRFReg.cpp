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
\brief Base class for all SIRF registration.

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFReg.h"
#include "SIRFRegMisc.h"
#include <_reg_localTransformation.h>
#include <nifti1_io.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <boost/filesystem.hpp>

using namespace std;

void SIRFReg::check_parameters()
{
    // If anything is missing
    if (_parameter_filename == "") {
        throw std::runtime_error("Parameter file has not been set.");}
    if (!_floating_image_sptr && _floating_image_filename == "") {
        throw std::runtime_error("Floating image has not been set."); }
    if (!_reference_image_sptr && _reference_image_filename == "") {
        throw std::runtime_error("Reference image has not been set."); }
}

void SIRFReg::save_warped_image(const string filename) const
{
    SIRFRegMisc::save_nifti_image(_warped_image_sptr,filename);
}

void SIRFReg::save_displacement_field_image(const string filename, bool split_xyz, bool flip_for_stir)
{
    shared_ptr<nifti_image> im_to_save_sptr;

    // If the user wants to save it as it is
    if (!flip_for_stir)
        im_to_save_sptr = _displacement_field_image_sptr;

    // But if they want to output for stir, they need to flip the z-axis
    else {
        im_to_save_sptr = make_shared<nifti_image>();
        // Deep copy the displacement image
        SIRFRegMisc::copy_nifti_image(im_to_save_sptr,_displacement_field_image_sptr);
        // Flip the z-axis
        SIRFRegMisc::flip_multicomponent_image(im_to_save_sptr,2);
    }

    // Whether anything has been flipped or not, save the result...

    // If the user wants it saved as multicomponent image
    if (!split_xyz)
        SIRFRegMisc::save_nifti_image(im_to_save_sptr,filename);

    // If the user wants the multicomponent image split into 3 separate images.
    else
        SIRFRegMisc::save_split_multicomponent_nifti_image(im_to_save_sptr,filename);
}

void SIRFReg::get_disp_from_cpp(shared_ptr<nifti_image> &cpp_sptr)
{
    // Calculate deformation field image
    nifti_image *def_ptr;
    def_ptr = nifti_copy_nim_info(_reference_image_sptr.get());
    def_ptr->dim[0]=def_ptr->ndim=5;
    def_ptr->dim[1]=def_ptr->nx=_reference_image_sptr->nx;
    def_ptr->dim[2]=def_ptr->ny=_reference_image_sptr->ny;
    def_ptr->dim[3]=def_ptr->nz=_reference_image_sptr->nz;
    def_ptr->dim[4]=def_ptr->nt=1;def_ptr->pixdim[4]=def_ptr->dt=1.0;
    if(_reference_image_sptr->nz>1) def_ptr->dim[5]=def_ptr->nu=3;
    else def_ptr->dim[5]=def_ptr->nu=2;
    def_ptr->pixdim[5]=def_ptr->du=1.0;
    def_ptr->dim[6]=def_ptr->nv=1;def_ptr->pixdim[6]=def_ptr->dv=1.0;
    def_ptr->dim[7]=def_ptr->nw=1;def_ptr->pixdim[7]=def_ptr->dw=1.0;
    def_ptr->nvox=def_ptr->nx*def_ptr->ny*def_ptr->nz*def_ptr->nt*def_ptr->nu;
    def_ptr->datatype = cpp_sptr->datatype;
    def_ptr->nbyper = cpp_sptr->nbyper;
    def_ptr->data = (void *)calloc(def_ptr->nvox, def_ptr->nbyper);
    reg_spline_getDeformationField(cpp_sptr.get(),
                                   _reference_image_sptr.get(),
                                   def_ptr,
                                   NULL,
                                   false, //composition
                                   true // bspline
                                   );

    shared_ptr<nifti_image> def_sptr = make_shared<nifti_image>(*def_ptr);

    // Get the disp field from the def field
    _displacement_field_image_sptr = make_shared<nifti_image>();
    SIRFRegMisc::copy_nifti_image(_displacement_field_image_sptr,def_sptr);

    reg_getDisplacementFromDeformation(_displacement_field_image_sptr.get());
}
