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

#include "SIRFRegActivityCorrect.h"
#include "SIRFRegMisc.h"
#include "_reg_tools.h"

using namespace std;

SIRFRegActivityCorrect::SIRFRegActivityCorrect()
{
    // Set all the default values
    _initial_activity   = -1.;
    _half_life          = -1.;
    _start              = -1.;
    _stop               = -1.;
}

void SIRFRegActivityCorrect::check_parameters()
{
    // If anything is missing
    if (_initial_activity == -1.) {
        throw std::runtime_error("Initial activity has not been set.");}
    if (_half_life == -1.) {
        throw std::runtime_error("Half lifehas not been set.");}
    if (_start == -1.) {
        throw std::runtime_error("Start time has not been set.");}
    if (_stop == -1.) {
        throw std::runtime_error("Stop time has not been set.");}
    if (!_input_image_sptr && _input_image_filename == "") {
        throw std::runtime_error("Floating image has not been set."); }
}

void SIRFRegActivityCorrect::update()
{
    cout << "\n\nCorrecting activity...\n\n";

    this->check_parameters();

    // Open images if necessary, correct if not
    if (!_input_image_sptr) {
        SIRFRegMisc::open_nifti_image(_input_image_sptr,_input_image_filename); }
    else {
        reg_checkAndCorrectDimension(_input_image_sptr.get()); }

    // Create the output as a copy of the input
    SIRFRegMisc::copy_nifti_image(_output_image_sptr,_input_image_sptr);

    // Calculate lambda
    float lambda = log(2) / _half_life;

    // ac_cumu = A_0/lambda * (exp(-lambda*start) - exp(-lambda*stop))
    float ac_cumu = ( _initial_activity / lambda )  *  (  exp(-lambda*_start)  -  exp(-lambda*_stop)  );

    // Do the division (last arg: 0=add,1=sub,2=mul,3=div)
    reg_tools_addSubMulDivValue(_input_image_sptr.get(),_output_image_sptr.get(),ac_cumu,3);
}

void SIRFRegActivityCorrect::save_output(const string filename) const
{
    SIRFRegMisc::save_nifti_image(_output_image_sptr,filename);
}
