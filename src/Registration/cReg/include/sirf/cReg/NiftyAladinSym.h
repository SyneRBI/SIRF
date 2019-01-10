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

#ifndef _SIRFREGNIFTYALADINSYM_H_
#define _SIRFREGNIFTYALADINSYM_H_

#include "sirf/cReg/NiftyRegistration.h"

template<class dataType> class reg_aladin_sym;

namespace sirf {

/// Forward declarations
template<class dataType> class AffineTransformation;

/// Wrapper around NiftyReg's aladin class for rigid and affine transformations
template<class dataType> class NiftyAladinSym : public NiftyRegistration<dataType>
{
public:

    /// Process
    void process();

    /// Get forwards transformation matrix
    const std::shared_ptr<const AffineTransformation<dataType> > get_transformation_matrix_forward() const { return _TM_forward_sptr; }

    /// Get inverse transformation matrix
    const std::shared_ptr<const AffineTransformation<dataType> > get_transformation_matrix_inverse() const { return _TM_inverse_sptr; }

protected:

    /// Parse parameter file
    virtual void parse_parameter_file();

    /// Set extra parameters.
    void set_parameters();

    /// Register object
    std::shared_ptr<reg_aladin_sym<dataType> > _registration_sptr;

    /// Forwards transformation matrix
    std::shared_ptr<AffineTransformation<dataType> > _TM_forward_sptr;
    /// Inverse transformation matrix
    std::shared_ptr<AffineTransformation<dataType> > _TM_inverse_sptr;
};
}

#endif
