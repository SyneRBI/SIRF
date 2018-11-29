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
#include "SIRFRegParser.h"
#include "NiftiImageData3DTensor.h"
#include <_reg_f3d_sym.h>
#include <_reg_base.h>

using namespace sirf;

template<class T>
void SIRFRegNiftyF3dSym<T>::process()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // Create the registration object
    _registration_sptr = std::shared_ptr<reg_f3d_sym<T> >(new reg_f3d_sym<T>(_reference_time_point, _floating_time_point));
    _registration_sptr->SetFloatingImage(_floating_image.get_raw_nifti_sptr().get());
    _registration_sptr->SetReferenceImage(_reference_image.get_raw_nifti_sptr().get());

    // If there is an initial transformation matrix, set it
    if (_use_initial_transformation) {
        mat44 init_tm = _initial_transformation.get_as_mat44();
        _registration_sptr->SetAffineTransformation(&init_tm);
    }

    // Set masks
    if (_reference_mask.is_initialised())
        _registration_sptr->SetReferenceMask(_reference_mask.get_raw_nifti_sptr().get());
    if (_floating_mask.is_initialised())
        _registration_sptr->SetFloatingMask(_floating_mask.get_raw_nifti_sptr().get());

    // Parse parameter file
    this->parse_parameter_file();

    // Set any extra parameters
    this->set_parameters();

    std::cout << "\n\nStarting registration...\n\n";

    // Run
    _registration_sptr->Run();

    // Get the warped image
    _warped_image = NiftiImageData3D(**_registration_sptr->GetWarpedImage());

    // For some reason, dt & pixdim[4] are sometimes set to 1
    if (_floating_image.get_raw_nifti_sptr()->dt < 1.e-7F &&
            _reference_image.get_raw_nifti_sptr()->dt < 1.e-7F)
        _warped_image.get_raw_nifti_sptr()->pixdim[4] = _warped_image.get_raw_nifti_sptr()->dt = 0.F;

    // Get the CPP images
    NiftiImageData3DTensor cpp_forward(*_registration_sptr->GetControlPointPositionImage());
    NiftiImageData3DTensor cpp_inverse(*_registration_sptr->GetBackwardControlPointPositionImage());

    // Get deformation fields from cpp
    _def_image_forward.create_from_cpp(cpp_forward, _reference_image);
    _def_image_inverse.create_from_cpp(cpp_inverse, _reference_image);

    // Get the displacement fields from the def
    _disp_image_forward.create_from_def(_def_image_forward);
    _disp_image_inverse.create_from_def(_def_image_inverse);

    std::cout << "\n\nRegistration finished!\n\n";
}

template<class T>
void SIRFRegNiftyF3dSym<T>::check_parameters()
{
    SIRFReg::check_parameters();

    // If anything is missing
    if (_floating_time_point == -1) {
        throw std::runtime_error("Floating time point has not been set."); }
    if (_reference_time_point == -1) {
        throw std::runtime_error("Reference time point has not been set."); }
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
    parser.add_key      ( "SetKLDWeight",                   &reg_f3d_sym<T>::SetKLDWeight                   );
    parser.add_key      ( "SetLinearEnergyWeight",          &reg_f3d_sym<T>::SetLinearEnergyWeight          );
    parser.add_key      ( "SetLNCCKernelType",              &reg_f3d_sym<T>::SetLNCCKernelType              );
    parser.add_key      ( "SetLNCCWeight",                  &reg_f3d_sym<T>::SetLNCCWeight                  );
    parser.add_key      ( "SetNMIWeight",                   &reg_f3d_sym<T>::SetNMIWeight                   );
    parser.add_key      ( "SetPerturbationNumber",          &reg_f3d_sym<T>::SetPerturbationNumber          );
    parser.add_key      ( "SetSSDWeight",                   &reg_f3d_sym<T>::SetSSDWeight                   );
    parser.parse();
}
template<class T>
void SIRFRegNiftyF3dSym<T>::set_parameters()
{
    for (size_t i=0; i<_extra_params.size(); i+=3) {

        std::string par  = _extra_params[ i ];
        std::string arg1 = _extra_params[i+1];
        std::string arg2 = _extra_params[i+2];

        // 1 argument
        if      (strcmp(par.c_str(),"SetBendingEnergyWeight")       == 0) _registration_sptr->SetBendingEnergyWeight        (stof(arg1));
        else if (strcmp(par.c_str(),"SetCompositionStepNumber")     == 0) _registration_sptr->SetCompositionStepNumber      (stoi(arg1));
        else if (strcmp(par.c_str(),"SetFloatingSmoothingSigma")    == 0) _registration_sptr->SetFloatingSmoothingSigma     (stof(arg1));
        else if (strcmp(par.c_str(),"SetGradientSmoothingSigma")    == 0) _registration_sptr->SetGradientSmoothingSigma     (stof(arg1));
        else if (strcmp(par.c_str(),"SetInverseConsistencyWeight")  == 0) _registration_sptr->SetInverseConsistencyWeight   (stof(arg1));
        else if (strcmp(par.c_str(),"SetJacobianLogWeight")         == 0) _registration_sptr->SetJacobianLogWeight          (stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceSmoothingSigma")   == 0) _registration_sptr->SetReferenceSmoothingSigma    (stof(arg1));
        else if (strcmp(par.c_str(),"SetWarpedPaddingValue")        == 0) _registration_sptr->SetWarpedPaddingValue         (stof(arg1));
        else if (strcmp(par.c_str(),"SetLinearEnergyWeight")        == 0) _registration_sptr->SetLinearEnergyWeight         (stof(arg1));
        else if (strcmp(par.c_str(),"SetLNCCKernelType")            == 0) _registration_sptr->SetLNCCKernelType             (stoi(arg1));
        else if (strcmp(par.c_str(),"SetLevelNumber")               == 0) _registration_sptr->SetLevelNumber                (unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetLevelToPerform")            == 0) _registration_sptr->SetLevelToPerform             (unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetMaximalIterationNumber")    == 0) _registration_sptr->SetMaximalIterationNumber     (unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetPerturbationNumber")        == 0) _registration_sptr->SetPerturbationNumber         (unsigned(stoi(arg1)));

        // 2 arguments
        else if (strcmp(par.c_str(),"SetSSDWeight")                 == 0) _registration_sptr->SetSSDWeight                  (stoi(arg1), stoi(arg2));
        else if (strcmp(par.c_str(),"SetLNCCWeight")                == 0) _registration_sptr->SetLNCCWeight                 (stoi(arg1), stod(arg1));
        else if (strcmp(par.c_str(),"SetNMIWeight")                 == 0) _registration_sptr->SetNMIWeight                  (stoi(arg1), stod(arg2));
        else if (strcmp(par.c_str(),"SetKLDWeight")                 == 0) _registration_sptr->SetKLDWeight                  (stoi(arg1), unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetFloatingThresholdUp")       == 0) _registration_sptr->SetFloatingThresholdUp        (unsigned(stoi(arg1)), stof(arg2));
        else if (strcmp(par.c_str(),"SetFloatingThresholdLow")      == 0) _registration_sptr->SetFloatingThresholdLow       (unsigned(stoi(arg1)), stoi(arg2));
        else if (strcmp(par.c_str(),"SetReferenceThresholdLow")     == 0) _registration_sptr->SetReferenceThresholdLow      (unsigned(stoi(arg1)), stof(arg2));
        else if (strcmp(par.c_str(),"SetReferenceThresholdUp")      == 0) _registration_sptr->SetReferenceThresholdUp       (unsigned(stoi(arg1)), stof(arg2));
        else if (strcmp(par.c_str(),"SetSpacing")                   == 0) _registration_sptr->SetSpacing                    (unsigned(stoi(arg1)), stof(arg2));

        else
            throw std::runtime_error("\nUnknown argument: " + par);
    }
}

namespace sirf {
//template class SIRFRegNiftyF3dSym<double>;
template class SIRFRegNiftyF3dSym<float>;
}
