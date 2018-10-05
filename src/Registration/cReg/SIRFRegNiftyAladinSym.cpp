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
\brief NiftyReg's aladin class for rigid and affine registrations.

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegNiftyAladinSym.h"
#include "SIRFRegMisc.h"
#include "SIRFRegParser.h"
#include <_reg_localTrans.h>
#include <_reg_tools.h>

using namespace std;
using namespace sirf;

template<class T>
void SIRFRegNiftyAladinSym<T>::update()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // Create the registration object
    _registration_sptr = make_shared<reg_aladin_sym<T> >();
    _registration_sptr->SetInputReference(_reference_image.get_raw_nifti_sptr().get());
    _registration_sptr->SetInputFloating(_floating_image.get_raw_nifti_sptr().get());

    // Parse parameter file
    this->parse_parameter_file();

    // Set any extra parameters
    this->set_parameters();

    cout << "\n\nStarting registration...\n\n";

    // Run
    _registration_sptr->Run();

    // Get the output
    _warped_image = NiftiImage3D(*_registration_sptr->GetFinalWarpedImage());

    // For some reason, dt & pixdim[4] are sometimes set to 1
    if (_floating_image.get_raw_nifti_sptr()->dt < 1.e-7F &&
            _reference_image.get_raw_nifti_sptr()->dt < 1.e-7F)
        _warped_image.get_raw_nifti_sptr()->pixdim[4] = _warped_image.get_raw_nifti_sptr()->dt = 0.F;

    // Get the forward and backward transformation matrices
    _TM_fwrd = *_registration_sptr->GetTransformationMatrix();
    _TM_back = nifti_mat44_inverse(_TM_fwrd.get_raw_mat44());

    cout << "\nPrinting forwards tranformation matrix:\n";
    _TM_fwrd.print();
    cout << "\nPrinting backwards tranformation matrix:\n";
    _TM_back.print();

    // affine->def->disp
    _def_image_fwrd.create_from_3D_image(_reference_image);
    _def_image_back.create_from_3D_image(_reference_image);

    mat44 temp_fwrd(_TM_fwrd.get_raw_mat44());
    mat44 temp_back(_TM_back.get_raw_mat44());
    reg_affine_getDeformationField(&temp_fwrd, _def_image_fwrd.get_raw_nifti_sptr().get());
    reg_affine_getDeformationField(&temp_back, _def_image_back.get_raw_nifti_sptr().get());

    // Get the displacement fields from the def
    _disp_image_fwrd.create_from_def(_def_image_fwrd);
    _disp_image_back.create_from_def(_def_image_back);

    cout << "\n\nRegistration finished!\n\n";
}

template<class T>
void SIRFRegNiftyAladinSym<T>::parse_parameter_file()
{
    SIRFRegParser<reg_aladin_sym<T> > parser;

    parser.set_object   ( _registration_sptr  );
    parser.set_filename ( _parameter_filename );
    parser.add_key      ( "SetAlignCentre",                     &reg_aladin<T>::SetAlignCentre                      );
    parser.add_key      ( "SetBlockPercentage",                 &reg_aladin<T>::SetBlockPercentage                  );
    parser.add_key      ( "SetFloatingSigma",                   &reg_aladin<T>::SetFloatingSigma                    );
    parser.add_key      ( "SetInlierLts",                       &reg_aladin<T>::SetInlierLts                        );
    parser.add_key      ( "SetInputTransform",                  &reg_aladin<T>::SetInputTransform                   );
    parser.add_key      ( "SetInterpolation",                   &reg_aladin<T>::SetInterpolation                    );
    parser.add_key      ( "SetInterpolationToCubic",            &reg_aladin<T>::SetInterpolationToCubic             );
    parser.add_key      ( "SetInterpolationToNearestNeighbor",  &reg_aladin<T>::SetInterpolationToNearestNeighbor   );
    parser.add_key      ( "SetInterpolationToTrilinear",        &reg_aladin<T>::SetInterpolationToTrilinear         );
    parser.add_key      ( "SetLevelsToPerform",                 &reg_aladin<T>::SetLevelsToPerform                  );
    parser.add_key      ( "SetMaxIterations",                   &reg_aladin<T>::SetMaxIterations                    );
    parser.add_key      ( "SetNumberOfLevels",                  &reg_aladin<T>::SetNumberOfLevels                   );
    parser.add_key      ( "SetPerformAffine",                   &reg_aladin<T>::SetPerformAffine                    );
    parser.add_key      ( "SetPerformRigid",                    &reg_aladin<T>::SetPerformRigid                     );
    parser.add_key      ( "SetReferenceSigma",                  &reg_aladin<T>::SetReferenceSigma                   );
    parser.add_key      ( "SetAlignCentreGravity",              &reg_aladin_sym<T>::SetAlignCentreGravity           );
    parser.add_key      ( "SetBlockStepSize",                   &reg_aladin_sym<T>::SetBlockStepSize                );
    parser.add_key      ( "SetFloatingLowerThreshold",          &reg_aladin_sym<T>::SetFloatingLowerThreshold       );
    parser.add_key      ( "SetFloatingUpperThreshold",          &reg_aladin_sym<T>::SetFloatingUpperThreshold       );
    parser.add_key      ( "SetReferenceLowerThreshold",         &reg_aladin_sym<T>::SetReferenceLowerThreshold      );
    parser.add_key      ( "SetReferenceUpperThreshold",         &reg_aladin_sym<T>::SetReferenceUpperThreshold      );
    parser.add_key      ( "SetVerbose",                         &reg_aladin_sym<T>::SetVerbose                      );
    parser.add_key      ( "SetWarpedPaddingValue",              &reg_aladin_sym<T>::SetWarpedPaddingValue           );
    parser.add_key      ( "setCaptureRangeVox",                 &reg_aladin_sym<T>::setCaptureRangeVox              );
    parser.add_key      ( "setGpuIdx",                          &reg_aladin_sym<T>::setGpuIdx                       );
    parser.add_key      ( "setPlatformCode",                    &reg_aladin_sym<T>::setPlatformCode                 );
    parser.parse();
}

template<class T>
void SIRFRegNiftyAladinSym<T>::set_parameters()
{
    for (size_t i=0; i<_extra_params.size(); i+=3) {

        string par  = _extra_params[ i ];
        string arg1 = _extra_params[i+1];
        string arg2 = _extra_params[i+2];

        // Void
        if      (strcmp(par.c_str(),"SetInterpolationToCubic")           == 0) _registration_sptr->SetInterpolationToCubic();
        else if (strcmp(par.c_str(),"SetInterpolationToNearestNeighbor") == 0) _registration_sptr->SetInterpolationToNearestNeighbor();
        else if (strcmp(par.c_str(),"SetInterpolationToTrilinear")       == 0) _registration_sptr->SetInterpolationToTrilinear();

        // String
        else if (strcmp(par.c_str(),"SetAlignCentre")                    == 0) _registration_sptr->SetAlignCentre(arg1.c_str());
        else if (strcmp(par.c_str(),"SetInputTransform")                 == 0) _registration_sptr->SetInputTransform(arg1.c_str());
        else if (strcmp(par.c_str(),"SetPerformAffine")                  == 0) _registration_sptr->SetPerformAffine(arg1.c_str());
        else if (strcmp(par.c_str(),"SetPerformRigid")                   == 0) _registration_sptr->SetPerformRigid(arg1.c_str());
        else if (strcmp(par.c_str(),"SetAlignCentreGravity")             == 0) _registration_sptr->SetAlignCentreGravity(arg1.c_str());
        else if (strcmp(par.c_str(),"SetVerbose")                        == 0) _registration_sptr->SetVerbose(arg1.c_str());

        // Int
        else if (strcmp(par.c_str(),"SetBlockPercentage")                == 0) _registration_sptr->SetBlockPercentage(stoi(arg1));
        else if (strcmp(par.c_str(),"SetInterpolation")                  == 0) _registration_sptr->SetInterpolation(stoi(arg1));
        else if (strcmp(par.c_str(),"SetBlockStepSize")                  == 0) _registration_sptr->SetBlockStepSize(stoi(arg1));
        else if (strcmp(par.c_str(),"setCaptureRangeVox")                == 0) _registration_sptr->setCaptureRangeVox(stoi(arg1));
        else if (strcmp(par.c_str(),"setPlatformCode")                   == 0) _registration_sptr->setPlatformCode(stoi(arg1));

        // Unsigned
        else if (strcmp(par.c_str(),"SetLevelsToPerform")                == 0) _registration_sptr->SetLevelsToPerform(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetMaxIterations")                  == 0) _registration_sptr->SetMaxIterations(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetNumberOfLevels")                 == 0) _registration_sptr->SetNumberOfLevels(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"setGpuIdx")                         == 0) _registration_sptr->setGpuIdx(unsigned(stoi(arg1)));

        // Float
        else if (strcmp(par.c_str(),"SetFloatingSigma")                  == 0) _registration_sptr->SetFloatingSigma(stof(arg1));
        else if (strcmp(par.c_str(),"SetInlierLts")                      == 0) _registration_sptr->SetInlierLts(stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceSigma")                 == 0) _registration_sptr->SetReferenceSigma(stof(arg1));
        else if (strcmp(par.c_str(),"SetFloatingLowerThreshold")         == 0) _registration_sptr->SetFloatingLowerThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetFloatingUpperThreshold")         == 0) _registration_sptr->SetFloatingUpperThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceLowerThreshold")        == 0) _registration_sptr->SetReferenceLowerThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceUpperThreshold")        == 0) _registration_sptr->SetReferenceUpperThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetWarpedPaddingValue")             == 0) _registration_sptr->SetWarpedPaddingValue(stof(arg1));

        else
            throw runtime_error("\nUnknown argument: " + par);
    }
}

namespace sirf {
// Put the instantiations of the template class at the END of the file!
template class SIRFRegNiftyAladinSym<float>;
template class SIRFRegNiftyAladinSym<double>;
}
