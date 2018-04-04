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

#ifndef _SIRFREGNIFTYALADIN_H_
#define _SIRFREGNIFTYALADIN_H_


#include "SIRFReg.h"
template<class T> class reg_aladin;

/// Wrapper around NiftyReg's aladin class for rigid and affine transformations
template<class T> class SIRFRegNiftyAladin : public SIRFReg
{
public:

    /// Update
    void update();

    /// Get transformation matrix
    mat44 *get_transformation_matrix() const { return _transformation_matrix_sptr.get(); }

    /// Get inverse transformation matrix
    mat44 *get_inverse_transformation_matrix() const { return _transformation_matrix_inverse_sptr.get(); }

    /// Save transformation matrix to file
    void save_transformation_matrix(const std::string filename) const;

    /// Save inverse transformation matrix to file
    void save_inverse_transformation_matrix(const std::string filename) const;

protected:

    /// Parse parameter file
    virtual void parse_parameter_file();

    /// Set up CPP image
    void set_up_CPP(std::shared_ptr<nifti_image> &cpp_sptr);

    /// Register object
    std::shared_ptr<reg_aladin<T> > _registration_sptr;

    /// Transformation matrix
    std::shared_ptr<mat44>          _transformation_matrix_sptr;
    /// Inverse transformation matrix
    std::shared_ptr<mat44>          _transformation_matrix_inverse_sptr;
};

#endif
