/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2019 University College London

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

#pragma once

#include "sirf/cReg/NiftyRegistration.h"

template<class dataType> class reg_f3d_sym;

namespace sirf {

/// Forward declarations
template<class dataType> class AffineTransformation;

/*!
\ingroup Registration
\brief NiftyReg's f3d class for non-rigid registrations.

User can set an initial affine transformation if desired. 

In theory, multiple time points can be used, but thus far has only been tested for
t == 1 for both reference and floating images.

\author Richard Brown
\author CCP PETMR
*/
template<class dataType> class NiftyF3dSym : public NiftyRegistration<dataType>
{
public:

    /// Constructor
    NiftyF3dSym()
    {
        _floating_time_point  = 1;
        _reference_time_point = 1;
    }

    /// Process
    void process();

    /// Set floating time point
    void set_floating_time_point(const int floating_time_point) { _floating_time_point = floating_time_point; }

    /// Set reference time point
    void set_reference_time_point(const int reference_time_point) { _reference_time_point = reference_time_point; }

    /// Set initial affine transformation
    void set_initial_affine_transformation(const std::shared_ptr<const AffineTransformation<float> > mat) { _initial_transformation_sptr = mat; }

protected:

    /// Check parameters
    virtual void check_parameters() const;

    /// Parse parameter file
    virtual void parse_parameter_file();

    /// Set extra parameters.
    void set_parameters();

    /// Registration object
    std::shared_ptr<reg_f3d_sym<dataType> > _registration_sptr;

    /// Floating time point
    int _floating_time_point;
    /// Reference time point
    int _reference_time_point;
    /// Transformation matrix
    std::shared_ptr<const AffineTransformation<float> > _initial_transformation_sptr;
};
}
