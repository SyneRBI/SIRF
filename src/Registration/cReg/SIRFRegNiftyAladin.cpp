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

#include "SIRFRegNiftyAladin.h"
#include "SIRFRegMisc.h"
#include "SIRFRegParser.h"
#include <_reg_aladin.h>
#include <_reg_tools.h>
#include <_reg_localTrans.h>

using namespace std;

template<class T>
void SIRFRegNiftyAladin<T>::update()
{
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
    _registration_sptr = make_shared<reg_aladin<T> >();
    _registration_sptr->SetInputReference(_reference_image_sptr.get());
    _registration_sptr->SetInputFloating(_floating_image_sptr.get());

    // Parse parameter file
    this->parse_parameter_file();

    cout << "\n\nStarting registration...\n\n";

    // Run
    _registration_sptr->Run();

    // Get the output
    _warped_image_sptr = std::make_shared<nifti_image>(*_registration_sptr->GetFinalWarpedImage());

    // Get the forward and backward transformation matrices
    _TM_fwrd_sptr = std::make_shared<mat44>(*_registration_sptr->GetTransformationMatrix());
    _TM_back_sptr = std::make_shared<mat44>(nifti_mat44_inverse(*_TM_fwrd_sptr.get()));

    cout << "\nPrinting forwards tranformation matrix:\n";
    SIRFRegMisc::print_mat44(_TM_fwrd_sptr.get());
    cout << "\nPrinting backwards tranformation matrix:\n";
    SIRFRegMisc::print_mat44(_TM_back_sptr.get());

    // Need to prepare a deformation field image for getting the forwards and backwards disp images
    shared_ptr<nifti_image> def_sptr;
    SIRFRegMisc::create_def_or_disp_image(def_sptr,_reference_image_sptr);

    // forward: affine->def->disp
    reg_affine_getDeformationField(_TM_fwrd_sptr.get(), def_sptr.get());
    _disp_image_fwrd_sptr = make_shared<nifti_image>();
    SIRFRegMisc::copy_nifti_image(_disp_image_fwrd_sptr,def_sptr);
    reg_getDisplacementFromDeformation(_disp_image_fwrd_sptr.get());

    // backward: affine->def->disp
    reg_affine_getDeformationField(_TM_back_sptr.get(), def_sptr.get());
    _disp_image_back_sptr = make_shared<nifti_image>();
    SIRFRegMisc::copy_nifti_image(_disp_image_back_sptr,def_sptr);
    reg_getDisplacementFromDeformation(_disp_image_back_sptr.get());

    cout << "\n\nRegistration finished!\n\n";
}

template<class T>
void SIRFRegNiftyAladin<T>::parse_parameter_file()
{
    SIRFRegParser<reg_aladin<T> > parser;
    parser.set_object   ( _registration_sptr  );
    parser.set_filename ( _parameter_filename );
    parser.add_key      ( "SetAlignCentre",                     &reg_aladin<T>::SetAlignCentre                      );
    parser.add_key      ( "SetAlignCentreGravity",              &reg_aladin<T>::SetAlignCentreGravity               );
    parser.add_key      ( "SetBlockPercentage",                 &reg_aladin<T>::SetBlockPercentage                  );
    parser.add_key      ( "SetBlockStepSize",                   &reg_aladin<T>::SetBlockStepSize                    );
    parser.add_key      ( "SetFloatingLowerThreshold",          &reg_aladin<T>::SetFloatingLowerThreshold           );
    parser.add_key      ( "SetFloatingSigma",                   &reg_aladin<T>::SetFloatingSigma                    );
    parser.add_key      ( "SetFloatingUpperThreshold",          &reg_aladin<T>::SetFloatingUpperThreshold           );
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
    parser.add_key      ( "SetReferenceLowerThreshold",         &reg_aladin<T>::SetReferenceLowerThreshold          );
    parser.add_key      ( "SetReferenceSigma",                  &reg_aladin<T>::SetReferenceSigma                   );
    parser.add_key      ( "SetReferenceUpperThreshold",         &reg_aladin<T>::SetReferenceUpperThreshold          );
    parser.add_key      ( "SetVerbose",                         &reg_aladin<T>::SetVerbose                          );
    parser.add_key      ( "SetWarpedPaddingValue",              &reg_aladin<T>::SetWarpedPaddingValue               );
    parser.add_key      ( "setGpuIdx",                          &reg_aladin<T>::setGpuIdx                           );
    parser.add_key      ( "setPlatformCode",                    &reg_aladin<T>::setPlatformCode                     );

    parser.parse();
}

template<class T>
void SIRFRegNiftyAladin<T>::save_transformation_matrix_fwrd(const std::string &filename) const
{
    save_transformation_matrix(_TM_fwrd_sptr,filename);
}

template<class T>
void SIRFRegNiftyAladin<T>::save_transformation_matrix_back(const std::string &filename) const
{
    save_transformation_matrix(_TM_back_sptr,filename);
}

template<class T>
void SIRFRegNiftyAladin<T>::
save_transformation_matrix(const std::shared_ptr<mat44> &TM_sptr, const std::string &filename) const
{
    // Check that the matrix exists
    if (!TM_sptr)
        throw std::runtime_error("Transformation matrix not available. Have you run the registration?");

    // Check that input isn't blank
    if (filename == "")
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");

    cout << "\nSaving inverse transformation matrix to file (" << filename << ")..." << flush;

    reg_tool_WriteAffineFile(TM_sptr.get(), filename.c_str());

    cout << "Done.\n";
}

// Put the instantiations of the template class at the END of the file!
template class SIRFRegNiftyAladin<float>;
template class SIRFRegNiftyAladin<double>;
