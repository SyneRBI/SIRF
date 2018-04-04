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
\brief Resampling class based on nifty resample

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegNiftyResample.h"
#include "SIRFRegMisc.h"
#include <_reg_resampling.h>
#include <_reg_globalTransformation.h>
#include <_reg_tools.h>

using namespace std;

void SIRFRegNiftyResample::update()
{
    cout << "\n\nStarting resampling...\n\n";

    // Check that all the required information has been entered
    check_parameters();

    // Open images if necessary, correct if not
    if (!_reference_image_sptr) {
        SIRFRegMisc::open_nifti_image(_reference_image_sptr,_reference_image_filename); }
    else {
        reg_checkAndCorrectDimension(_reference_image_sptr.get()); }

    if (!_floating_image_sptr) {
        SIRFRegMisc::open_nifti_image(_floating_image_sptr,_floating_image_filename); }
    else {
        reg_checkAndCorrectDimension(_floating_image_sptr.get()); }

    // Set up transformation matrix
    mat44 transformation_matrix;
    set_up_transformation_matrix(transformation_matrix);

    cout << "\n\nthe transformation matrix is:\n";
    SIRFRegMisc::print_mat44(&transformation_matrix);

    cout << "\n\nConverting affine transformation to deformation field...\n\n";

    // Initialise the deformation field image
    shared_ptr<nifti_image> deformation_field_image_sptr;
    this->set_up_deformation_field_image(deformation_field_image_sptr, transformation_matrix);

    cout << "\n\nSuccessfully converted affine transformation to deformation field.\n\n";

    // Setup output image
    set_up_output_image();

    reg_resampleSourceImage(_reference_image_sptr.get(),
                                _floating_image_sptr.get(),
                                _output_image_sptr.get(),
                                deformation_field_image_sptr.get(),
                                NULL,
                                _interpolation_type,
                                0);

    cout << "\n\nResampling finished!\n\n";
}

void SIRFRegNiftyResample::check_parameters()
{
    // If anything is missing
    if (!_reference_image_sptr && _reference_image_filename == "") {
        throw std::runtime_error("Reference image has not been set."); }
    if (!_floating_image_sptr && _floating_image_filename == "") {
        throw std::runtime_error("Floating image has not been set."); }
    if (_transformation_matrices.size() == 0) {
        throw std::runtime_error("Transformation matrix/matrices not set."); }
    if (_interpolation_type == NOTSET) {
        throw std::runtime_error("Interpolation type has not been set."); }
}

void SIRFRegNiftyResample::save_resampled_image(const string filename) const
{
    SIRFRegMisc::save_nifti_image(_output_image_sptr,filename);
}

void SIRFRegNiftyResample::set_up_transformation_matrix(mat44 &matrix)
{
    // Start as identity
    for (int i=0;i<4;i++) {
        for (int j=0; j<4; j++) {
            if (i==j) matrix.m[i][j] = 1.;
            else      matrix.m[i][j] = 0.;
        }
    }

    // Loop over the transformation matrices
    for (int i=0; i<_transformation_matrices.size(); i++) {
        matrix = SIRFRegMisc::multiply_mat44(matrix,*_transformation_matrices[i].get());
    }

    cout << "\n\nHere's the result of the matrix multiplications:\n";
    SIRFRegMisc::print_mat44(&matrix);
    cout << "\n\n";
}

void SIRFRegNiftyResample::set_up_deformation_field_image(shared_ptr<nifti_image> &deformation_field_image_sptr, mat44 matrix)
{
    // Copy the info from the reference image
    deformation_field_image_sptr = shared_ptr<nifti_image>(nifti_copy_nim_info(_reference_image_sptr.get()));

    // Set up the rest of the info
    deformation_field_image_sptr->dim[0]=deformation_field_image_sptr->ndim=5;
    deformation_field_image_sptr->dim[1]=deformation_field_image_sptr->nx=_reference_image_sptr->nx;
    deformation_field_image_sptr->dim[2]=deformation_field_image_sptr->ny=_reference_image_sptr->ny;
    deformation_field_image_sptr->dim[3]=deformation_field_image_sptr->nz=_reference_image_sptr->nz;
    deformation_field_image_sptr->dim[4]=deformation_field_image_sptr->nt=1;deformation_field_image_sptr->pixdim[4]=deformation_field_image_sptr->dt=1.0;
    if(_reference_image_sptr->nz>1) deformation_field_image_sptr->dim[5]=deformation_field_image_sptr->nu=3;
    else deformation_field_image_sptr->dim[5]=deformation_field_image_sptr->nu=2;
    deformation_field_image_sptr->pixdim[5]=deformation_field_image_sptr->du=1.0;
    deformation_field_image_sptr->dim[6]=deformation_field_image_sptr->nv=1;deformation_field_image_sptr->pixdim[6]=deformation_field_image_sptr->dv=1.0;
    deformation_field_image_sptr->dim[7]=deformation_field_image_sptr->nw=1;deformation_field_image_sptr->pixdim[7]=deformation_field_image_sptr->dw=1.0;
    deformation_field_image_sptr->nvox=deformation_field_image_sptr->nx*deformation_field_image_sptr->ny*deformation_field_image_sptr->nz*deformation_field_image_sptr->nt*deformation_field_image_sptr->nu;
    deformation_field_image_sptr->datatype = NIFTI_TYPE_FLOAT32; // need to add in an if/else NIFTI_TYPE_FLOAT64 if you want to template the class for doubles
    deformation_field_image_sptr->nbyper = sizeof(float);
    deformation_field_image_sptr->data = (void *)calloc(deformation_field_image_sptr->nvox, deformation_field_image_sptr->nbyper);

    reg_affine_positionField(&matrix,
                             _reference_image_sptr.get(),
                             deformation_field_image_sptr.get());
}

void SIRFRegNiftyResample::set_up_output_image()
{
    _output_image_sptr = shared_ptr<nifti_image>(nifti_copy_nim_info(_reference_image_sptr.get()));
    _output_image_sptr->dim[0]=_output_image_sptr->ndim=_floating_image_sptr->dim[0];
    _output_image_sptr->dim[4]=_output_image_sptr->nt=_floating_image_sptr->dim[4];
    _output_image_sptr->cal_min=_floating_image_sptr->cal_min;
    _output_image_sptr->cal_max=_floating_image_sptr->cal_max;
    _output_image_sptr->scl_slope=_floating_image_sptr->scl_slope;
    _output_image_sptr->scl_inter=_floating_image_sptr->scl_inter;
    _output_image_sptr->datatype = _floating_image_sptr->datatype;
    _output_image_sptr->nbyper = _floating_image_sptr->nbyper;
    _output_image_sptr->nvox = _output_image_sptr->dim[1] * _output_image_sptr->dim[2] * _output_image_sptr->dim[3] * _output_image_sptr->dim[4];
    _output_image_sptr->data = (void *)calloc(_output_image_sptr->nvox, _output_image_sptr->nbyper);
}
