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
#include <_reg_globalTrans.h>
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
    SIRFRegMisc::create_def_or_disp_image(deformation_field_image_sptr,_reference_image_sptr);
    reg_affine_getDeformationField(&transformation_matrix,deformation_field_image_sptr.get());

    SIRFRegMisc::save_nifti_image(deformation_field_image_sptr,"/Users/rich/Desktop/temp/def");

    cout << "\n\nSuccessfully converted affine transformation to deformation field.\n\n";

    // Setup output image
    _output_image_sptr = std::make_shared<nifti_image>();
    SIRFRegMisc::copy_nifti_image(_output_image_sptr,_floating_image_sptr);

    reg_resampleImage(_reference_image_sptr.get(),
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
