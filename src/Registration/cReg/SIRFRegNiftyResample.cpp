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
#include "NiftiImageData3DTensor.h"
#include "NiftiImageData3DDeformation.h"
#include "NiftiImageData3DDisplacement.h"
#include "SIRFRegMat44.h"
#include <_reg_resampling.h>
#include <_reg_globalTrans.h>
#include <_reg_tools.h>
#include <memory>

using namespace sirf;

void SIRFRegNiftyResample::add_transformation_affine(const SIRFRegMat44 &affine)
{
    _transformations.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegMat44(affine.deep_copy())));
}

void SIRFRegNiftyResample::add_transformation_disp(const NiftiImageData3DDisplacement &disp)
{
    _transformations.push_back(std::shared_ptr<SIRFRegTransformation>(new NiftiImageData3DDisplacement(disp.deep_copy())));
}

void SIRFRegNiftyResample::add_transformation_def(const NiftiImageData3DDeformation &def)
{
    _transformations.push_back(std::shared_ptr<SIRFRegTransformation>(new NiftiImageData3DDeformation(def.deep_copy())));
}

void SIRFRegNiftyResample::process()
{
    std::cout << "\n\nStarting resampling...\n\n";

    // Check that all the required information has been entered
    check_parameters();

    // Compose single transformation from multiple
    NiftiImageData3DDeformation transformation =
            NiftiImageData3DDeformation::compose_single_deformation(_transformations, _reference_image);

    // Setup output image
    set_up_output_image();

    reg_resampleImage(_floating_image.get_raw_nifti_sptr().get(),
                      _output_image.get_raw_nifti_sptr().get(),
                      transformation.get_as_deformation_field(_reference_image).get_raw_nifti_sptr().get(),
                      NULL,
                      _interpolation_type,
                      0);

    is_finished = true;

    std::cout << "\n\nResampling finished!\n\n";
}

const NiftiImageData3D &SIRFRegNiftyResample::get_output() const
{
    if (!is_finished)
        throw std::runtime_error("Cannot return output image until resampling has been performed.");
    return _output_image;
}

void SIRFRegNiftyResample::check_parameters()
{
    // If anything is missing
    if (!_reference_image.is_initialised()) {
        throw std::runtime_error("Reference image has not been set."); }
    if (!_floating_image.is_initialised()) {
        throw std::runtime_error("Floating image has not been set."); }

    if (_transformations.size() == 0)
        throw std::runtime_error("Transformation(s) not set.");
}

void SIRFRegNiftyResample::set_up_output_image()
{
    _output_image = _reference_image;

    const nifti_image *floating_ptr = _floating_image.get_raw_nifti_sptr().get();
    nifti_image       *output_ptr   = _output_image.get_raw_nifti_sptr().get();

    output_ptr->cal_min                   = floating_ptr->cal_min;
    output_ptr->cal_max                   = floating_ptr->cal_max;
    output_ptr->scl_slope                 = floating_ptr->scl_slope;
    output_ptr->scl_inter                 = floating_ptr->scl_inter;
    output_ptr->datatype = floating_ptr->datatype;
    output_ptr->intent_code=floating_ptr->intent_code;
    memset(output_ptr->intent_name, 0, 16);
    strcpy(output_ptr->intent_name,floating_ptr->intent_name);
    output_ptr->intent_p1 = floating_ptr->intent_p1;
    output_ptr->intent_p2 = floating_ptr->intent_p2;
    output_ptr->datatype  = floating_ptr->datatype;
    output_ptr->nbyper    = floating_ptr->nbyper;
    output_ptr->nvox = unsigned(output_ptr->dim[1] * output_ptr->dim[2] * output_ptr->dim[3] * output_ptr->dim[4] * output_ptr->dim[5]);
    output_ptr->data = static_cast<void *>(calloc(output_ptr->nvox, unsigned(output_ptr->nbyper)));
}
