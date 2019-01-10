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

#include "sirf/cReg/NiftyAladinSym.h"
#include "sirf/cReg/Parser.h"
#include "sirf/cReg/AffineTransformation.h"
#include "sirf/cReg/NiftiImageData3D.h"
#include "sirf/cReg/NiftiImageData3DDeformation.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"
#include <_reg_aladin_sym.h>

using namespace sirf;

template<class dataType>
void NiftyAladinSym<dataType>::process()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // Convert the input images from ImageData to NiftiImageData3D
    this->set_up_inputs();

    // Annoyingly NiftyReg doesn't mark ref and floating images as const, so need to copy (could do a naughty cast, but not going to do that!)
    NiftiImageData3D<dataType> ref = *this->_reference_image_nifti_sptr;
    NiftiImageData3D<dataType> flo = *this->_floating_image_nifti_sptr;

    // Create the registration object
    _registration_sptr = std::make_shared<reg_aladin_sym<dataType> >();
    _registration_sptr->SetInputReference(ref.get_raw_nifti_sptr().get());
    _registration_sptr->SetInputFloating(flo.get_raw_nifti_sptr().get());

    // Set masks (if present). Again, need to copy to get rid of const
    if (this->_reference_mask_nifti_sptr && this->_reference_mask_nifti_sptr->is_initialised()) {
        NiftiImageData3D<dataType> ref_mask = *this->_reference_mask_nifti_sptr;
        _registration_sptr->SetInputMask(ref_mask.get_raw_nifti_sptr().get());
    }
    if (this->_floating_mask_nifti_sptr && this->_floating_mask_nifti_sptr->is_initialised()) {
        NiftiImageData3D<dataType> flo_mask = *this->_floating_mask_nifti_sptr;
        _registration_sptr->SetInputMask(flo_mask.get_raw_nifti_sptr().get());
    }

    // Parse parameter file
    this->parse_parameter_file();

    // Set any extra parameters
    this->set_parameters();

    std::cout << "\n\nStarting registration...\n\n";

    // Run
    _registration_sptr->Run();

    // Get the output
    this->_warped_image_nifti_sptr = std::make_shared<NiftiImageData3D<dataType> >(*_registration_sptr->GetFinalWarpedImage());

    // For some reason, dt & pixdim[4] are sometimes set to 1
    if (this->_floating_image_nifti_sptr->get_raw_nifti_sptr()->dt < 1.e-7F &&
            this->_reference_image_nifti_sptr->get_raw_nifti_sptr()->dt < 1.e-7F)
        this->_warped_image_nifti_sptr->get_raw_nifti_sptr()->pixdim[4] = this->_warped_image_nifti_sptr->get_raw_nifti_sptr()->dt = 0.F;

    // Get the forward and inverse transformation matrices
    this->_TM_forward_sptr = std::make_shared<AffineTransformation<dataType> >(*this->_registration_sptr->GetTransformationMatrix());
    this->_TM_inverse_sptr = std::make_shared<AffineTransformation<dataType> >(nifti_mat44_inverse(*_registration_sptr->GetTransformationMatrix()));

    std::cout << "\nPrinting forwards tranformation matrix:\n";
    _TM_forward_sptr->print();
    std::cout << "\nPrinting inverse tranformation matrix:\n";
    std::dynamic_pointer_cast<const AffineTransformation<dataType> >(_TM_inverse_sptr)->print();

    // Get as deformation and displacement
    NiftiImageData3DDeformation<dataType> def_fwd = _TM_forward_sptr->get_as_deformation_field(ref);
    NiftiImageData3DDeformation<dataType> def_inv = _TM_inverse_sptr->get_as_deformation_field(ref);
    this->_disp_image_forward_sptr = std::make_shared<NiftiImageData3DDisplacement<dataType> >(def_fwd);
    this->_disp_image_inverse_sptr = std::make_shared<NiftiImageData3DDisplacement<dataType> >(def_inv);

    // The output should be a clone of the reference image, with data filled in from the nifti image
    this->_warped_image_sptr = this->_reference_image_sptr->clone();
    this->_warped_image_sptr->fill(*this->_warped_image_nifti_sptr);

    std::cout << "\n\nRegistration finished!\n\n";
}

template<class dataType>
void NiftyAladinSym<dataType>::parse_parameter_file()
{
    Parser<reg_aladin_sym<dataType> > parser;

    parser.set_object   (    _registration_sptr     );
    parser.set_filename ( this->_parameter_filename );
    parser.add_key      ( "SetAlignCentre",                     &reg_aladin<dataType>::SetAlignCentre                      );
    parser.add_key      ( "SetBlockPercentage",                 &reg_aladin<dataType>::SetBlockPercentage                  );
    parser.add_key      ( "SetFloatingSigma",                   &reg_aladin<dataType>::SetFloatingSigma                    );
    parser.add_key      ( "SetInlierLts",                       &reg_aladin<dataType>::SetInlierLts                        );
    parser.add_key      ( "SetInputTransform",                  &reg_aladin<dataType>::SetInputTransform                   );
    parser.add_key      ( "SetInterpolation",                   &reg_aladin<dataType>::SetInterpolation                    );
    parser.add_key      ( "SetInterpolationToCubic",            &reg_aladin<dataType>::SetInterpolationToCubic             );
    parser.add_key      ( "SetInterpolationToNearestNeighbor",  &reg_aladin<dataType>::SetInterpolationToNearestNeighbor   );
    parser.add_key      ( "SetInterpolationToTrilinear",        &reg_aladin<dataType>::SetInterpolationToTrilinear         );
    parser.add_key      ( "SetLevelsToPerform",                 &reg_aladin<dataType>::SetLevelsToPerform                  );
    parser.add_key      ( "SetMaxIterations",                   &reg_aladin<dataType>::SetMaxIterations                    );
    parser.add_key      ( "SetNumberOfLevels",                  &reg_aladin<dataType>::SetNumberOfLevels                   );
    parser.add_key      ( "SetPerformAffine",                   &reg_aladin<dataType>::SetPerformAffine                    );
    parser.add_key      ( "SetPerformRigid",                    &reg_aladin<dataType>::SetPerformRigid                     );
    parser.add_key      ( "SetReferenceSigma",                  &reg_aladin<dataType>::SetReferenceSigma                   );
    parser.add_key      ( "SetAlignCentreGravity",              &reg_aladin_sym<dataType>::SetAlignCentreGravity           );
    parser.add_key      ( "SetBlockStepSize",                   &reg_aladin_sym<dataType>::SetBlockStepSize                );
    parser.add_key      ( "SetFloatingLowerThreshold",          &reg_aladin_sym<dataType>::SetFloatingLowerThreshold       );
    parser.add_key      ( "SetFloatingUpperThreshold",          &reg_aladin_sym<dataType>::SetFloatingUpperThreshold       );
    parser.add_key      ( "SetReferenceLowerThreshold",         &reg_aladin_sym<dataType>::SetReferenceLowerThreshold      );
    parser.add_key      ( "SetReferenceUpperThreshold",         &reg_aladin_sym<dataType>::SetReferenceUpperThreshold      );
    parser.add_key      ( "SetVerbose",                         &reg_aladin_sym<dataType>::SetVerbose                      );
    parser.add_key      ( "SetWarpedPaddingValue",              &reg_aladin_sym<dataType>::SetWarpedPaddingValue           );
    parser.add_key      ( "setCaptureRangeVox",                 &reg_aladin_sym<dataType>::setCaptureRangeVox              );
    parser.add_key      ( "setGpuIdx",                          &reg_aladin_sym<dataType>::setGpuIdx                       );
    parser.add_key      ( "setPlatformCode",                    &reg_aladin_sym<dataType>::setPlatformCode                 );
    parser.parse();
}

template<class dataType>
void NiftyAladinSym<dataType>::set_parameters()
{
    for (size_t i=0; i<this->_extra_params.size(); i+=3) {

        std::string par  = this->_extra_params[ i ];
        std::string arg1 = this->_extra_params[i+1];
        // std::string arg2 = _extra_params[i+2]; No aladin methods need 2 args (but f3d does)

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
            throw std::runtime_error("\nUnknown argument: " + par);
    }
}

namespace sirf {
template class NiftyAladinSym<float>;
}
