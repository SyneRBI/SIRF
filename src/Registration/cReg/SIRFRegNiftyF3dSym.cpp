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
    if (!_reference_image_sptr)
        SIRFRegMisc::open_nifti_image(_reference_image_sptr,_reference_image_filename);
    else
        reg_checkAndCorrectDimension(_reference_image_sptr.get());

    if (!_floating_image_sptr)
        SIRFRegMisc::open_nifti_image(_floating_image_sptr,_floating_image_filename);
    else
        reg_checkAndCorrectDimension(_floating_image_sptr.get());

    // Create the registration object
    _registration_sptr = std::shared_ptr<reg_f3d_sym<T> >(new reg_f3d_sym<T>(_reference_time_point, _floating_time_point));
    _registration_sptr->SetFloatingImage(_floating_image_sptr.get());
    _registration_sptr->SetReferenceImage(_reference_image_sptr.get());

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
    _warped_image_sptr = std::shared_ptr<nifti_image>(*_registration_sptr->GetWarpedImage());

    // Get the CPP images
    std::shared_ptr<nifti_image> cpp_fwrd_sptr = std::shared_ptr<nifti_image>(_registration_sptr->GetControlPointPositionImage());
    std::shared_ptr<nifti_image> cpp_back_sptr = std::shared_ptr<nifti_image>(_registration_sptr->GetBackwardControlPointPositionImage());

    // Need to correct the CPP images (otherwise nv=0 and you can't read with matlab)
    reg_checkAndCorrectDimension(cpp_fwrd_sptr.get());
    reg_checkAndCorrectDimension(cpp_back_sptr.get());

    // Get deformation fields from cpp
    SIRFRegMisc::get_def_from_cpp(_def_image_fwrd_sptr,cpp_fwrd_sptr, _reference_image_sptr);
    SIRFRegMisc::get_def_from_cpp(_def_image_back_sptr,cpp_back_sptr, _reference_image_sptr);

    // Get the displacement fields from the def
    SIRFRegMisc::get_disp_from_def(_disp_image_fwrd_sptr,_def_image_fwrd_sptr);
    SIRFRegMisc::get_disp_from_def(_disp_image_back_sptr,_def_image_back_sptr);

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
