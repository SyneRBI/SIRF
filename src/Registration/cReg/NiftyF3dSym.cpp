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
\brief NiftyReg's f3d class for non-rigid registrations.

\author Richard Brown
\author SyneRBI
*/

#include "sirf/Reg/NiftyF3dSym.h"
#include "sirf/Reg/Parser.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include <_reg_f3d_sym.h>

using namespace sirf;

template<class dataType>
void NiftyF3dSym<dataType>::process()
{
    // Check the paramters that are NOT set via the parameter file have been set.
    this->check_parameters();

    // Convert the input images from ImageData to NiftiImageData3D
    this->set_up_inputs();

    // Annoyingly NiftyReg doesn't mark ref and floating images as const, so need to copy (could do a naughty cast, but not going to do that!)
    NiftiImageData3D<dataType> ref = *this->_reference_image_nifti_sptr;
    NiftiImageData3D<dataType> flo = *this->_floating_images_nifti.at(0);

    // Create the registration object
    if (_use_symmetric)
        _registration_sptr = std::make_shared<reg_f3d_sym<dataType> >(_reference_time_point, _floating_time_point);
    else
        _registration_sptr = std::make_shared<reg_f3d<dataType> >(_reference_time_point, _floating_time_point);

    // Set reference and floating images
    _registration_sptr->SetReferenceImage(ref.get_raw_nifti_sptr().get());
    _registration_sptr->SetFloatingImage(flo.get_raw_nifti_sptr().get());

    // By default, use a padding value of 0
    _registration_sptr->SetWarpedPaddingValue(0.f);

    // If there is an initial transformation matrix, set it
    if (_initial_transformation_sptr) {
        mat44 init_tm = _initial_transformation_sptr->get_as_mat44();
        _registration_sptr->SetAffineTransformation(&init_tm);
    }

    // Set masks (if present). Again, need to copy to get rid of const
    NiftiImageData3D<dataType> ref_mask, flo_mask;
    if (this->_reference_mask_nifti_sptr && this->_reference_mask_nifti_sptr->is_initialised()) {
        ref_mask = *this->_reference_mask_nifti_sptr;
        _registration_sptr->SetReferenceMask(ref_mask.get_raw_nifti_sptr().get());
    }
    if (this->_floating_mask_nifti_sptr && this->_floating_mask_nifti_sptr->is_initialised()) {
        flo_mask = *this->_floating_mask_sptr;
        _registration_sptr->SetFloatingMask(flo_mask.get_raw_nifti_sptr().get());
    }

    // Parse parameter file
    this->parse_parameter_file();

    // Set any extra parameters
    this->set_parameters();

    std::cout << "\n\nStarting registration...\n\n";

    // Run
    _registration_sptr->Run();

    // Get the warped image
    nifti_image **warped_im = _registration_sptr->GetWarpedImage();
    this->_warped_images_nifti.at(0) = std::make_shared<NiftiImageData3D<dataType> >(*warped_im[0]);
    // Free the images created
    if(warped_im[0]!=NULL)
        nifti_image_free(warped_im[0]);
    warped_im[0]=NULL;
    if(warped_im[1]!=NULL)
        nifti_image_free(warped_im[1]);
    warped_im[1]=NULL;
    free(warped_im);
    warped_im=NULL;

    // For some reason, dt & pixdim[4] are sometimes set to 1
    if (this->_floating_images_nifti.at(0)->get_raw_nifti_sptr()->dt < 1.e-7F &&
            this->_reference_image_nifti_sptr->get_raw_nifti_sptr()->dt < 1.e-7F)
        this->_warped_images_nifti.at(0)->get_raw_nifti_sptr()->pixdim[4] = this->_warped_images_nifti.at(0)->get_raw_nifti_sptr()->dt = 0.F;

    // Get the CPP images
    nifti_image * cpp_fwd_ptr = _registration_sptr->GetControlPointPositionImage();
    NiftiImageData3DTensor<dataType> cpp_forward(*cpp_fwd_ptr);
    nifti_image_free(cpp_fwd_ptr);

    // Get deformation fields from cpp
    std::shared_ptr<NiftiImageData3DDeformation<dataType> > def_fwd_sptr = std::make_shared<NiftiImageData3DDeformation<dataType> >();
    def_fwd_sptr->create_from_cpp(cpp_forward, ref);
    this->_def_fwd_images.at(0) = def_fwd_sptr;

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
NiftyF3dSym<dataType>::
get_deformation_field_inverse_sptr(const unsigned idx) const
{
    if (idx>0)
        throw std::runtime_error("NiftyF3dSym::get_displacement_field_inverse_sptr: idx out of range");

    std::shared_ptr<NiftiImageData3DDeformation<dataType> > def_inv_sptr;
    std::shared_ptr<Transformation<dataType> > trans_fwd = this->_def_fwd_images.at(0);
    const NiftiImageData3DDeformation<dataType> &def_fwd =
            *std::dynamic_pointer_cast<NiftiImageData3DDeformation<dataType> >(trans_fwd);

    // Get inverse deformation.
    // NiftyReg can only do inverse for 3D images.
    if (def_fwd.get_raw_nifti_sptr()->nu == 3)
        def_inv_sptr = def_fwd.get_inverse(this->_floating_images_nifti.at(0));
    else {
#ifdef SIRF_VTK
        // if not 3d but VTK is present, use that.
        def_inv_sptr = def_fwd.get_inverse(nullptr, true);
#else
        // else, throw an error
        throw std::runtime_error("F3D: Inversion not currently implemented for 2D images with NiftyReg. "
                                 "This means that you won't be able to get the inverse displacement or "
                                 "deformation fields. To have these, use a 3D image or install VTK");
#endif
    }
    return def_inv_sptr;
}

template<class dataType>
void NiftyF3dSym<dataType>::print_all_wrapped_methods()
{
    std::cout << "SetCompositionStepNumber(int)\n"
                 "SetInverseConsistencyWeight(dataType)\n"
                 "SetJacobianLogWeight(dataType)\n"
                 "SetLinearEnergyWeight(dataType)\n"
                 "SetWarpedPaddingValue(dataType)\n"
                 "SetBendingEnergyWeight(dataType)\n"
                 "SetFloatingSmoothingSigma(dataType)\n"
                 "SetGradientSmoothingSigma(dataType)\n"
                 "SetReferenceSmoothingSigma(dataType)\n"
                 "SetLNCCKernelType(int)\n"
                 "SetLevelNumber(unsigned)\n"
                 "SetLevelToPerform(unsigned)\n"
                 "SetMaximalIterationNumber(unsigned)\n"
                 "SetPerturbationNumber(unsigned)\n"
                 "SetSSDWeight(int,int)\n"
                 "SetLNCCWeight(int,double)\n"
                 "SetNMIWeight(int,double)\n"
                 "SetKLDWeight(int,unsigned)\n"
                 "SetFloatingThresholdUp(unsigned,dataType)\n"
                 "SetFloatingThresholdLow(unsigned,dataType)\n"
                 "SetReferenceThresholdUp(unsigned,dataType)\n"
                 "SetReferenceThresholdLow(unsigned,dataType)\n"
                 "SetSpacing(unsigned,dataType)\n";
}

template<class dataType>
void NiftyF3dSym<dataType>::check_parameters() const
{
    Registration<dataType>::check_parameters();

    // If anything is missing
    if (_floating_time_point == -1) {
        throw std::runtime_error("Floating time point has not been set."); }
    if (_reference_time_point == -1) {
        throw std::runtime_error("Reference time point has not been set."); }
}

template<class dataType>
void NiftyF3dSym<dataType>::parse_parameter_file()
{
    if (this->_parameter_filename.empty())
        return;

    Parser<reg_f3d<dataType> > parser;
    parser.set_object   (    _registration_sptr     );
    parser.set_filename ( this->_parameter_filename );

    parser.add_key("SetCompositionStepNumber",&reg_f3d<dataType>::SetCompositionStepNumber);
    parser.add_key("SetInverseConsistencyWeight",&reg_f3d<dataType>::SetInverseConsistencyWeight);
    parser.add_key("SetJacobianLogWeight",&reg_f3d<dataType>::SetJacobianLogWeight);
    parser.add_key("SetLinearEnergyWeight",&reg_f3d<dataType>::SetLinearEnergyWeight);
    parser.add_key("SetWarpedPaddingValue",&reg_f3d<dataType>::SetWarpedPaddingValue);
    parser.add_key("SetBendingEnergyWeight",&reg_f3d<dataType>::SetBendingEnergyWeight);
    parser.add_key("SetFloatingSmoothingSigma",&reg_f3d<dataType>::SetFloatingSmoothingSigma);
    parser.add_key("SetGradientSmoothingSigma",&reg_f3d<dataType>::SetGradientSmoothingSigma);
    parser.add_key("SetReferenceSmoothingSigma",&reg_f3d<dataType>::SetReferenceSmoothingSigma);
    parser.add_key("SetLNCCKernelType",&reg_f3d<dataType>::SetLNCCKernelType);
    parser.add_key("SetLevelNumber",&reg_f3d<dataType>::SetLevelNumber);
    parser.add_key("SetLevelToPerform",&reg_f3d<dataType>::SetLevelToPerform);
    parser.add_key("SetMaximalIterationNumber",&reg_f3d<dataType>::SetMaximalIterationNumber);
    parser.add_key("SetPerturbationNumber",&reg_f3d<dataType>::SetPerturbationNumber);
    parser.add_key("SetSSDWeight",&reg_f3d<dataType>::SetSSDWeight);
    parser.add_key("SetLNCCWeight",&reg_f3d<dataType>::SetLNCCWeight);
    parser.add_key("SetNMIWeight",&reg_f3d<dataType>::SetNMIWeight);
    parser.add_key("SetKLDWeight",&reg_f3d<dataType>::SetKLDWeight);
    parser.add_key("SetFloatingThresholdUp",&reg_f3d<dataType>::SetFloatingThresholdUp);
    parser.add_key("SetFloatingThresholdLow",&reg_f3d<dataType>::SetFloatingThresholdLow);
    parser.add_key("SetReferenceThresholdUp",&reg_f3d<dataType>::SetReferenceThresholdUp);
    parser.add_key("SetReferenceThresholdLow",&reg_f3d<dataType>::SetReferenceThresholdLow);
    parser.add_key("SetSpacing",&reg_f3d<dataType>::SetSpacing);

    parser.parse();
}
template<class dataType>
void NiftyF3dSym<dataType>::set_parameters()
{
    for (size_t i=0; i<this->_extra_params.size(); i+=3) {

        std::string par  = this->_extra_params[ i ];
        std::string arg1 = this->_extra_params[i+1];
        std::string arg2 = this->_extra_params[i+2];

        if      (strcmp(par.c_str(),"SetCompositionStepNumber")== 0) _registration_sptr->SetCompositionStepNumber(stoi(arg1));
        else if (strcmp(par.c_str(),"SetInverseConsistencyWeight")== 0) _registration_sptr->SetInverseConsistencyWeight(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetJacobianLogWeight")== 0) _registration_sptr->SetJacobianLogWeight(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetLinearEnergyWeight")== 0) _registration_sptr->SetLinearEnergyWeight(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetWarpedPaddingValue")== 0) _registration_sptr->SetWarpedPaddingValue(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetBendingEnergyWeight")== 0) _registration_sptr->SetBendingEnergyWeight(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetFloatingSmoothingSigma")== 0) _registration_sptr->SetFloatingSmoothingSigma(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetGradientSmoothingSigma")== 0) _registration_sptr->SetGradientSmoothingSigma(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetReferenceSmoothingSigma")== 0) _registration_sptr->SetReferenceSmoothingSigma(dataType(stod(arg1)));
        else if (strcmp(par.c_str(),"SetLNCCKernelType")== 0) _registration_sptr->SetLNCCKernelType(stoi(arg1));
        else if (strcmp(par.c_str(),"SetLevelNumber")== 0) _registration_sptr->SetLevelNumber(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetLevelToPerform")== 0) _registration_sptr->SetLevelToPerform(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetMaximalIterationNumber")== 0) _registration_sptr->SetMaximalIterationNumber(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetPerturbationNumber")== 0) _registration_sptr->SetPerturbationNumber(unsigned(stoi(arg1)));
        else if (strcmp(par.c_str(),"SetSSDWeight")== 0) _registration_sptr->SetSSDWeight(stoi(arg1), stoi(arg2));
        else if (strcmp(par.c_str(),"SetLNCCWeight")== 0) _registration_sptr->SetLNCCWeight(stoi(arg1), stod(arg2));
        else if (strcmp(par.c_str(),"SetNMIWeight")== 0) _registration_sptr->SetNMIWeight(stoi(arg1), stod(arg2));
        else if (strcmp(par.c_str(),"SetKLDWeight")== 0) _registration_sptr->SetKLDWeight(stoi(arg1), unsigned(stoi(arg2)));
        else if (strcmp(par.c_str(),"SetFloatingThresholdUp")== 0) _registration_sptr->SetFloatingThresholdUp(unsigned(stoi(arg1)), dataType(stod(arg2)));
        else if (strcmp(par.c_str(),"SetFloatingThresholdLow")== 0) _registration_sptr->SetFloatingThresholdLow(unsigned(stoi(arg1)), dataType(stod(arg2)));
        else if (strcmp(par.c_str(),"SetReferenceThresholdUp")== 0) _registration_sptr->SetReferenceThresholdUp(unsigned(stoi(arg1)), dataType(stod(arg2)));
        else if (strcmp(par.c_str(),"SetReferenceThresholdLow")== 0) _registration_sptr->SetReferenceThresholdLow(unsigned(stoi(arg1)), dataType(stod(arg2)));
        else if (strcmp(par.c_str(),"SetSpacing")== 0) _registration_sptr->SetSpacing(unsigned(stoi(arg1)), dataType(stod(arg2)));

        else
            throw std::runtime_error("\nUnknown argument: " + par);
    }
}

namespace sirf {
template class NiftyF3dSym<float>;
}
