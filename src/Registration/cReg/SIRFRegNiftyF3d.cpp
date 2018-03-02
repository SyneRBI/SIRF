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
\brief NiftyReg's f3d class for non-rigid registrations.

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegNiftyF3d.h"
#include "SIRFRegMisc.h"
#include "SIRFRegParser.h"
#include "_reg_f3d.h"

using namespace std;

template<class T>
void SIRFRegNiftyF3d<T>::update()
{
    // Try
    try {

        // Check the paramters that are NOT set via the parameter file have been set.
        this->check_parameters();

        // Open images if necessary, correct if not
        if (!_reference_image_sptr) {
            SIRFRegMisc::open_nifti_image(_reference_image_sptr,_reference_image_filename); }
        else {
            reg_checkAndCorrectDimension(_reference_image_sptr.get()); }

        if (!_floating_image_sptr) {
            SIRFRegMisc::open_nifti_image(_floating_image_sptr,_floating_image_filename); }
        else {
            reg_checkAndCorrectDimension(_floating_image_sptr.get()); }

        // Create the registration object
        _registration_sptr = std::shared_ptr<reg_f3d<T> >(new reg_f3d<T>(_reference_time_point, _floating_time_point));
        _registration_sptr->SetFloatingImage(_floating_image_sptr.get());
        _registration_sptr->SetReferenceImage(_reference_image_sptr.get());

        // Parse parameter file
        this->parse_parameter_file();

        cout << "\n\nStarting registration...\n\n";

        // Run
        _registration_sptr->Run_f3d();

        // Get the output
        _warped_image_sptr = std::shared_ptr<nifti_image>(*_registration_sptr->GetWarpedImage());

        cout << "\n\nRegistration finished!\n\n";

    // If there was an error, rethrow it.
    } catch (const std::runtime_error &error) {
        throw;
    }
}

template<class T>
void SIRFRegNiftyF3d<T>::check_parameters()
{
    SIRFReg::check_parameters();

    // If anything is missing
    if (_floating_time_point == -1) {
        throw std::runtime_error("Floating time point has not been set."); }
    if (_reference_time_point == -1) {
        throw std::runtime_error("Reference time point has not been set."); }
}

template<class T>
void SIRFRegNiftyF3d<T>::parse_parameter_file()
{
    SIRFRegParser<reg_f3d<T> > parser;
    parser.set_object   ( _registration_sptr  );
    parser.set_filename ( _parameter_filename );
    parser.add_key      ("SetAdditiveMC",                &reg_f3d<T>::SetAdditiveMC               );
    parser.add_key      ( "SetBendingEnergyWeight",      &reg_f3d<T>::SetBendingEnergyWeight      );
    parser.add_key      ( "SetCompositionStepNumber",    &reg_f3d<T>::SetCompositionStepNumber    );
    parser.add_key      ( "SetFloatingBinNumber",        &reg_f3d<T>::SetFloatingBinNumber        );
    parser.add_key      ( "SetFloatingSmoothingSigma",   &reg_f3d<T>::SetFloatingSmoothingSigma   );
    parser.add_key      ( "SetFloatingThresholdLow",     &reg_f3d<T>::SetFloatingThresholdLow     );
    parser.add_key      ( "SetFloatingThresholdUp",      &reg_f3d<T>::SetFloatingThresholdUp      );
    parser.add_key      ( "SetGradientSmoothingSigma",   &reg_f3d<T>::SetGradientSmoothingSigma   );
    parser.add_key      ( "SetInverseConsistencyWeight", &reg_f3d<T>::SetInverseConsistencyWeight );
    parser.add_key      ( "SetJacobianLogWeight",        &reg_f3d<T>::SetJacobianLogWeight        );
    parser.add_key      ( "SetL2NormDisplacementWeight", &reg_f3d<T>::SetL2NormDisplacementWeight );
    parser.add_key      ( "SetLevelNumber",              &reg_f3d<T>::SetLevelNumber              );
    parser.add_key      ( "SetLevelToPerform",           &reg_f3d<T>::SetLevelToPerform           );
    parser.add_key      ( "SetLinearEnergyWeights",      &reg_f3d<T>::SetLinearEnergyWeights      );
    parser.add_key      ( "SetMaximalIterationNumber",   &reg_f3d<T>::SetMaximalIterationNumber   );
    parser.add_key      ( "SetReferenceBinNumber",       &reg_f3d<T>::SetReferenceBinNumber       );
    parser.add_key      ( "SetReferenceSmoothingSigma",  &reg_f3d<T>::SetReferenceSmoothingSigma  );
    parser.add_key      ( "SetReferenceThresholdLow",    &reg_f3d<T>::SetReferenceThresholdLow    );
    parser.add_key      ( "SetReferenceThresholdUp",     &reg_f3d<T>::SetReferenceThresholdUp     );
    parser.add_key      ( "SetSpacing",                  &reg_f3d<T>::SetSpacing                  );
    parser.add_key      ( "SetWarpedPaddingValue",       &reg_f3d<T>::SetWarpedPaddingValue       );

    parser.parse();
}

template class SIRFRegNiftyF3d<double>;
template class SIRFRegNiftyF3d<float>;
