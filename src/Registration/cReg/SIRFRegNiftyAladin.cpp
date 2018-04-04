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
#include <_reg_localTransformation.h>

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

    // Get the transformation matrix and its inverse
    _transformation_matrix_sptr         = std::make_shared<mat44>(*_registration_sptr->GetTransformationMatrix());
    _transformation_matrix_inverse_sptr = std::make_shared<mat44>(nifti_mat44_inverse(*_transformation_matrix_sptr.get()));

    cout << "\nPrinting tranformation matrix:\n";
    SIRFRegMisc::print_mat44(_transformation_matrix_sptr.get());
    cout << "\nPrinting inverse tranformation matrix:\n";
    SIRFRegMisc::print_mat44(_transformation_matrix_inverse_sptr.get());

    // Convert transformation matrix to cpp image
    shared_ptr<nifti_image> control_point_position_image_sptr;
    set_up_CPP(control_point_position_image_sptr);

    // Get the disp field from the cpp image
    this->get_disp_from_cpp(control_point_position_image_sptr);

    vector<shared_ptr<nifti_image> > gf;
    gf.push_back(_warped_image_sptr);
    gf.push_back(_displacement_field_image_sptr);
    SIRFRegMisc::dump_nifti_info(gf);

    cout << "\n\nRegistration finished!\n\n";
}

template<class T>
void SIRFRegNiftyAladin<T>::parse_parameter_file()
{
    SIRFRegParser<reg_aladin<T> > parser;
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

    parser.parse();
}

template<class T>
void SIRFRegNiftyAladin<T>::save_transformation_matrix(const std::string filename) const
{
    // Check that the matrix exists
    if (!_transformation_matrix_sptr) {
        throw std::runtime_error("Transformation matrix not available. Have you run the registration?");
    }

    // Check that input isn't blank
    if (filename == "") {
        throw std::runtime_error("Error, cannot write transformation matrix to file because filename is blank");
    }

    cout << "\nSaving transformation matrix to file (" << filename << ")..." << flush;

    reg_tool_WriteAffineFile(_transformation_matrix_sptr.get(), filename.c_str());

    cout << "Done.\n";
}

template<class T>
void SIRFRegNiftyAladin<T>::save_inverse_transformation_matrix(const std::string filename) const
{
    // Check that the matrix exists
    if (!_transformation_matrix_inverse_sptr) {
        throw std::runtime_error("Inverse transformation matrix not available. Have you run the registration?");
    }

    // Check that input isn't blank
    if (filename == "") {
        throw std::runtime_error("Error, cannot write inverse transformation matrix to file because filename is blank");
    }

    cout << "\nSaving inverse transformation matrix to file (" << filename << ")..." << flush;

    reg_tool_WriteAffineFile(_transformation_matrix_inverse_sptr.get(), filename.c_str());

    cout << "Done.\n";
}

template<class T>
void SIRFRegNiftyAladin<T>::set_up_CPP(shared_ptr<nifti_image> &cpp_sptr)
{
    // Copy info from the reference image
    nifti_image *cpp_ptr = cpp_sptr.get();
    cpp_ptr = nifti_copy_nim_info(_warped_image_sptr.get());

    // Edit some of the information to make it a cpp image
    cpp_ptr->dim[0]   = cpp_ptr->ndim = 5;
    cpp_ptr->dim[5]   = cpp_ptr->nu   = 3;
    cpp_ptr->datatype = 16;
    cpp_ptr->nbyper   = 4;
    cpp_ptr->nvox    *= 3;

    // Allocate memory
    cpp_ptr->data=(void *)malloc(cpp_ptr->nvox*cpp_ptr->nbyper);

    // Convert affine transformation to cpp
    reg_bspline_initialiseControlPointGridWithAffine(_transformation_matrix_sptr.get(), cpp_ptr);

    // Need to correct the control point position image (otherwise nv=0 and you can't read with matlab)
    reg_checkAndCorrectDimension(cpp_ptr);

    // Copy output
    cpp_sptr = make_shared<nifti_image>(*cpp_ptr);
}

// Put the instantiations of the template class at the END of the file!
template class SIRFRegNiftyAladin<float>;
template class SIRFRegNiftyAladin<double>;
