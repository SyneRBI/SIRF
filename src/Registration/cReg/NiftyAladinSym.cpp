/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\author SyneRBI
*/

#include "sirf/Reg/NiftyAladinSym.h"
#include "sirf/Reg/Parser.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
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
    NiftiImageData3D<dataType> flo = *this->_floating_images_nifti.at(0);

    // Create the registration object
    _registration_sptr = std::make_shared<reg_aladin_sym<dataType> >();
    _registration_sptr->SetInputReference(ref.get_raw_nifti_sptr().get());
    _registration_sptr->SetInputFloating(flo.get_raw_nifti_sptr().get());

    // By default, use a padding value of 0
    _registration_sptr->SetWarpedPaddingValue(0.f);

    // Set masks (if present). Again, need to copy to get rid of const
    NiftiImageData3D<dataType> ref_mask, flo_mask;
    if (this->_reference_mask_nifti_sptr && this->_reference_mask_nifti_sptr->is_initialised()) {
        ref_mask = *this->_reference_mask_nifti_sptr;
        _registration_sptr->SetInputMask(ref_mask.get_raw_nifti_sptr().get());
    }
    if (this->_floating_mask_nifti_sptr && this->_floating_mask_nifti_sptr->is_initialised()) {
        flo_mask = *this->_floating_mask_nifti_sptr;
        _registration_sptr->SetInputFloatingMask(flo_mask.get_raw_nifti_sptr().get());
    }

    // Parse parameter file
    this->parse_parameter_file();

    // Set any extra parameters
    this->set_parameters();

    std::cout << "\n\nStarting registration...\n\n";

    // Run
    _registration_sptr->Run();

    // Get the output
    nifti_image *warped_im = _registration_sptr->GetFinalWarpedImage();
    this->_warped_images_nifti.at(0) = std::make_shared<NiftiImageData3D<dataType> >(*warped_im);
    nifti_image_free(warped_im);

    // For some reason, dt & pixdim[4] are sometimes set to 1
    if (this->_floating_images_nifti.at(0)->get_raw_nifti_sptr()->dt < 1.e-7F &&
            this->_reference_image_nifti_sptr->get_raw_nifti_sptr()->dt < 1.e-7F)
        this->_warped_images_nifti.at(0)->get_raw_nifti_sptr()->pixdim[4] = this->_warped_images_nifti.at(0)->get_raw_nifti_sptr()->dt = 0.F;

    // Get the forward and inverse transformation matrices
    this->_TM_forward_sptr = std::make_shared<AffineTransformation<float> >(*this->_registration_sptr->GetTransformationMatrix());
    this->_TM_inverse_sptr = std::make_shared<AffineTransformation<float> >(nifti_mat44_inverse(*_registration_sptr->GetTransformationMatrix()));

    std::cout << "\nPrinting forwards tranformation matrix:\n";
    _TM_forward_sptr->print();
    std::cout << "\nPrinting inverse tranformation matrix:\n";
    std::dynamic_pointer_cast<const AffineTransformation<float> >(_TM_inverse_sptr)->print();

    this->_def_fwd_images.at(0) = std::make_shared<NiftiImageData3DDeformation<dataType> >(_TM_forward_sptr->get_as_deformation_field(ref));

    // Output different dependent on whether image was set as object or via filename
    if (this->_reference_image_sptr) {
        // The output should be a clone of the reference image, with data filled in from the nifti image
        this->_warped_images.at(0) = this->_reference_image_sptr->clone();
        this->_warped_images.at(0)->fill(*this->_warped_images_nifti.at(0));
    }
    else
        this->_warped_images.at(0) = this->_warped_images_nifti.at(0);

    std::cout << "\n\nRegistration finished!\n\n";
}

template<class dataType>
const std::shared_ptr<const Transformation<dataType> >
NiftyAladinSym<dataType>::
get_deformation_field_inverse_sptr(const unsigned idx) const
{
    if (idx>0)
        throw std::runtime_error("NiftyAladinSym::get_deformation_field_inverse_sptr: idx out of range");

    return std::make_shared<NiftiImageData3DDeformation<dataType> >(
                _TM_inverse_sptr->get_as_deformation_field(*this->_floating_images_nifti.at(0)));
}

template<class dataType>
void NiftyAladinSym<dataType>::print_all_wrapped_methods()
{
    std::cout << "SetInterpolationToCubic()\n"
                 "SetInterpolationToNearestNeighbor()\n"
                 "SetInterpolationToTrilinear()\n"
                 "SetAlignCentre(bool)\n"
                 "SetInputTransform(str)\n"
                 "SetPerformAffine(bool)\n"
                 "SetPerformRigid(bool)\n"
                 "SetVerbose(bool)\n"
                 "SetBlockPercentage(int)\n"
                 "SetInterpolation(int)\n"
                 "SetBlockStepSize(int)\n"
                 "setCaptureRangeVox(int)\n"
                 "setPlatformCode(int)\n"
                 "SetLevelsToPerform(unsigned)\n"
                 "SetMaxIterations(unsigned)\n"
                 "SetNumberOfLevels(unsigned)\n"
                 "setGpuIdx(unsigned)\n"
                 "SetFloatingSigma(float)\n"
                 "SetInlierLts(float)\n"
                 "SetReferenceSigma(float)\n"
                 "SetFloatingLowerThreshold(float)\n"
                 "SetFloatingUpperThreshold(float)\n"
                 "SetReferenceLowerThreshold(float)\n"
                 "SetReferenceUpperThreshold(float)\n"
                 "SetWarpedPaddingValue(float)\n"
                 "SetAlignCentreMass(int)\n";
}

template<class dataType>
void NiftyAladinSym<dataType>::parse_parameter_file()
{
    if (this->_parameter_filename.empty())
        return;

    Parser<reg_aladin_sym<dataType> > parser;

    parser.set_object   (    _registration_sptr     );
    parser.set_filename ( this->_parameter_filename );

    parser.add_key("SetInterpolationToCubic",&reg_aladin_sym<dataType>::SetInterpolationToCubic);
    parser.add_key("SetInterpolationToNearestNeighbor",&reg_aladin_sym<dataType>::SetInterpolationToNearestNeighbor);
    parser.add_key("SetInterpolationToTrilinear",&reg_aladin_sym<dataType>::SetInterpolationToTrilinear);
    parser.add_key("SetAlignCentre",&reg_aladin_sym<dataType>::SetAlignCentre);
    parser.add_key("SetInputTransform",&reg_aladin_sym<dataType>::SetInputTransform);
    parser.add_key("SetPerformAffine",&reg_aladin_sym<dataType>::SetPerformAffine);
    parser.add_key("SetPerformRigid",&reg_aladin_sym<dataType>::SetPerformRigid);
    parser.add_key("SetVerbose",&reg_aladin_sym<dataType>::SetVerbose);
    parser.add_key("SetBlockPercentage",&reg_aladin_sym<dataType>::SetBlockPercentage);
    parser.add_key("SetInterpolation",&reg_aladin_sym<dataType>::SetInterpolation);
    parser.add_key("SetBlockStepSize",&reg_aladin_sym<dataType>::SetBlockStepSize);
    parser.add_key("setCaptureRangeVox",&reg_aladin_sym<dataType>::setCaptureRangeVox);
    parser.add_key("setPlatformCode",&reg_aladin_sym<dataType>::setPlatformCode);
    parser.add_key("SetLevelsToPerform",&reg_aladin_sym<dataType>::SetLevelsToPerform);
    parser.add_key("SetMaxIterations",&reg_aladin_sym<dataType>::SetMaxIterations);
    parser.add_key("SetNumberOfLevels",&reg_aladin_sym<dataType>::SetNumberOfLevels);
    parser.add_key("setGpuIdx",&reg_aladin_sym<dataType>::setGpuIdx);
    parser.add_key("SetFloatingSigma",&reg_aladin_sym<dataType>::SetFloatingSigma);
    parser.add_key("SetInlierLts",&reg_aladin_sym<dataType>::SetInlierLts);
    parser.add_key("SetReferenceSigma",&reg_aladin_sym<dataType>::SetReferenceSigma);
    parser.add_key("SetFloatingLowerThreshold",&reg_aladin_sym<dataType>::SetFloatingLowerThreshold);
    parser.add_key("SetFloatingUpperThreshold",&reg_aladin_sym<dataType>::SetFloatingUpperThreshold);
    parser.add_key("SetReferenceLowerThreshold",&reg_aladin_sym<dataType>::SetReferenceLowerThreshold);
    parser.add_key("SetReferenceUpperThreshold",&reg_aladin_sym<dataType>::SetReferenceUpperThreshold);
    parser.add_key("SetWarpedPaddingValue",&reg_aladin_sym<dataType>::SetWarpedPaddingValue);
    parser.add_key("SetAlignCentreMass",&reg_aladin_sym<dataType>::SetAlignCentreMass);

    parser.parse();
}

template<class dataType>
void NiftyAladinSym<dataType>::set_parameters()
{
    for (size_t i=0; i<this->_extra_params.size(); i+=3) {

        std::string par  = this->_extra_params[ i ];
        std::string arg1 = this->_extra_params[i+1];
        // std::string arg2 = this->_extra_params[i+2]; No aladin methods need 2 args (but f3d does)

        if      (strcmp(par.c_str(),"SetInterpolationToCubic")== 0) _registration_sptr->SetInterpolationToCubic();
        else if (strcmp(par.c_str(),"SetInterpolationToNearestNeighbor")== 0) _registration_sptr->SetInterpolationToNearestNeighbor();
        else if (strcmp(par.c_str(),"SetInterpolationToTrilinear")== 0) _registration_sptr->SetInterpolationToTrilinear();
        else if (strcmp(par.c_str(),"SetAlignCentre")== 0) _registration_sptr->SetAlignCentre(bool(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetInputTransform")== 0) _registration_sptr->SetInputTransform(arg1.c_str());
        else if (strcmp(par.c_str(),"SetPerformAffine")== 0) _registration_sptr->SetPerformAffine(bool(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetPerformRigid")== 0) _registration_sptr->SetPerformRigid(bool(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetVerbose")== 0) _registration_sptr->SetVerbose(bool(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetBlockPercentage")== 0) _registration_sptr->SetBlockPercentage(stoi(arg1));
        else if (strcmp(par.c_str(),"SetInterpolation")== 0) _registration_sptr->SetInterpolation(stoi(arg1));
        else if (strcmp(par.c_str(),"SetBlockStepSize")== 0) _registration_sptr->SetBlockStepSize(stoi(arg1));
        else if (strcmp(par.c_str(),"setCaptureRangeVox")== 0) _registration_sptr->setCaptureRangeVox(stoi(arg1));
        else if (strcmp(par.c_str(),"setPlatformCode")== 0) _registration_sptr->setPlatformCode(stoi(arg1));
        else if (strcmp(par.c_str(),"SetLevelsToPerform")== 0) _registration_sptr->SetLevelsToPerform(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetMaxIterations")== 0) _registration_sptr->SetMaxIterations(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetNumberOfLevels")== 0) _registration_sptr->SetNumberOfLevels(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"setGpuIdx")== 0) _registration_sptr->setGpuIdx(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetFloatingSigma")== 0) _registration_sptr->SetFloatingSigma(stof(arg1));
        else if (strcmp(par.c_str(),"SetInlierLts")== 0) _registration_sptr->SetInlierLts(stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceSigma")== 0) _registration_sptr->SetReferenceSigma(stof(arg1));
        else if (strcmp(par.c_str(),"SetFloatingLowerThreshold")== 0) _registration_sptr->SetFloatingLowerThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetFloatingUpperThreshold")== 0) _registration_sptr->SetFloatingUpperThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceLowerThreshold")== 0) _registration_sptr->SetReferenceLowerThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetReferenceUpperThreshold")== 0) _registration_sptr->SetReferenceUpperThreshold(stof(arg1));
        else if (strcmp(par.c_str(),"SetWarpedPaddingValue")== 0) _registration_sptr->SetWarpedPaddingValue(stof(arg1));
        else if (strcmp(par.c_str(),"SetAlignCentreMass")== 0) _registration_sptr->SetAlignCentreMass(stoi(arg1));
        else
            throw std::runtime_error("\nUnknown argument: " + par);
    }
}

namespace sirf {
template class NiftyAladinSym<float>;
}
