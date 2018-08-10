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

#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegMisc.h"
#include "SIRFRegParser.h"
#include "SIRFImageDataDeformation.h"
#include <_reg_f3d_sym.h>
#if NIFTYREG_VER_1_5
#include <_reg_base.h>
#endif

using namespace std;

template<class T>
void SIRFRegNiftyF3dSym<T>::update()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // Open images if necessary, correct if not
    reg_checkAndCorrectDimension(_reference_image.get_image_as_nifti().get());
    reg_checkAndCorrectDimension(_floating_image.get_image_as_nifti().get());

    // Create the registration object
    _registration_sptr = std::shared_ptr<reg_f3d_sym<T> >(new reg_f3d_sym<T>(_reference_time_point, _floating_time_point));
    _registration_sptr->SetFloatingImage(_floating_image.get_image_as_nifti().get());
    _registration_sptr->SetReferenceImage(_reference_image.get_image_as_nifti().get());

    // If an initial transformation matrix has been set via filename, open it
    if (_initial_transformation_filename != "")
        SIRFRegMisc::open_transformation_matrix(_initial_transformation_sptr,_initial_transformation_filename);

    // If there is an initial transformation matrix, set it
    if (_initial_transformation_sptr)
        _registration_sptr->SetAffineTransformation(_initial_transformation_sptr.get());

    // Parse parameter file
    this->parse_parameter_file();

    cout << "\n\nStarting registration...\n\n";

    // Run
#if NIFTYREG_VER_1_5
    _registration_sptr->Run();
#elif NIFTYREG_VER_1_3
    _registration_sptr->Run_f3d();
#endif

    // Get the warped image
    _warped_image = SIRFImageData(*_registration_sptr->GetWarpedImage());

    // Get the CPP images
    SIRFImageDataDeformation cpp_fwrd(_registration_sptr->GetControlPointPositionImage());
    SIRFImageDataDeformation cpp_back(_registration_sptr->GetBackwardControlPointPositionImage());

    // Get deformation fields from cpp
    SIRFRegMisc::get_def_from_cpp(_def_image_fwrd,cpp_fwrd.get_image_as_nifti(), _reference_image);
    SIRFRegMisc::get_def_from_cpp(_def_image_back,cpp_back.get_image_as_nifti(), _reference_image);

    // Get the displacement fields from the def
    _disp_image_fwrd = _def_image_fwrd;
    _disp_image_back = _def_image_back;
    SIRFRegMisc::convert_from_def_to_disp(_disp_image_fwrd);
    SIRFRegMisc::convert_from_def_to_disp(_disp_image_back);

    cout << "\n\nRegistration finished!\n\n";
}

template<class T>
void SIRFRegNiftyF3dSym<T>::check_parameters()
{
    SIRFReg::check_parameters();

    // If anything is missing
    if (_floating_time_point == -1) {
        throw runtime_error("Floating time point has not been set."); }
    if (_reference_time_point == -1) {
        throw runtime_error("Reference time point has not been set."); }
}

template<class T>
void SIRFRegNiftyF3dSym<T>::parse_parameter_file()
{
    SIRFRegParser<reg_f3d_sym<T> > parser;
    parser.set_object   ( _registration_sptr  );
    parser.set_filename ( _parameter_filename );

    parser.add_key      ( "SetBendingEnergyWeight",         &reg_f3d_sym<T>::SetBendingEnergyWeight         );
    parser.add_key      ( "SetCompositionStepNumber",       &reg_f3d_sym<T>::SetCompositionStepNumber       );
    parser.add_key      ( "SetFloatingSmoothingSigma",      &reg_f3d_sym<T>::SetFloatingSmoothingSigma      );
    parser.add_key      ( "SetFloatingThresholdLow",        &reg_f3d_sym<T>::SetFloatingThresholdLow        );
    parser.add_key      ( "SetFloatingThresholdUp",         &reg_f3d_sym<T>::SetFloatingThresholdUp         );
    parser.add_key      ( "SetGradientSmoothingSigma",      &reg_f3d_sym<T>::SetGradientSmoothingSigma      );
    parser.add_key      ( "SetInverseConsistencyWeight",    &reg_f3d_sym<T>::SetInverseConsistencyWeight    );
    parser.add_key      ( "SetJacobianLogWeight",           &reg_f3d_sym<T>::SetJacobianLogWeight           );
    parser.add_key      ( "SetLevelNumber",                 &reg_f3d_sym<T>::SetLevelNumber                 );
    parser.add_key      ( "SetLevelToPerform",              &reg_f3d_sym<T>::SetLevelToPerform              );
    parser.add_key      ( "SetMaximalIterationNumber",      &reg_f3d_sym<T>::SetMaximalIterationNumber      );
    parser.add_key      ( "SetReferenceSmoothingSigma",     &reg_f3d_sym<T>::SetReferenceSmoothingSigma     );
    parser.add_key      ( "SetReferenceThresholdLow",       &reg_f3d_sym<T>::SetReferenceThresholdLow       );
    parser.add_key      ( "SetReferenceThresholdUp",        &reg_f3d_sym<T>::SetReferenceThresholdUp        );
    parser.add_key      ( "SetSpacing",                     &reg_f3d_sym<T>::SetSpacing                     );
    parser.add_key      ( "SetWarpedPaddingValue",          &reg_f3d_sym<T>::SetWarpedPaddingValue          );

#if NIFTYREG_VER_1_5
    parser.add_key      ( "SetKLDWeight",                   &reg_f3d_sym<T>::SetKLDWeight                   );
    parser.add_key      ( "SetLinearEnergyWeight",          &reg_f3d_sym<T>::SetLinearEnergyWeight          );
    parser.add_key      ( "SetLNCCKernelType",              &reg_f3d_sym<T>::SetLNCCKernelType              );
    parser.add_key      ( "SetLNCCWeight",                  &reg_f3d_sym<T>::SetLNCCWeight                  );
    parser.add_key      ( "SetNMIWeight",                   &reg_f3d_sym<T>::SetNMIWeight                   );
    parser.add_key      ( "SetPerturbationNumber",          &reg_f3d_sym<T>::SetPerturbationNumber          );
    parser.add_key      ( "SetSSDWeight",                   &reg_f3d_sym<T>::SetSSDWeight                   );
#endif
    parser.parse();
}

//template class SIRFRegNiftyF3dSym<double>;
template class SIRFRegNiftyF3dSym<float>;
