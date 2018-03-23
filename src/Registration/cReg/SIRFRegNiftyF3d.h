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

#ifndef _SIRFREGNIFTYF3D_H_
#define _SIRFREGNIFTYF3D_H_

#include "SIRFReg.h"

template<class T> class reg_f3d;

/// Wrapper around NiftyReg's f3d class for non-rigid transformations
template<class T> class SIRFRegNiftyF3d : public SIRFReg
{
public:

    /// Constructor
    SIRFRegNiftyF3d()
    {
        _floating_time_point  = -1;
        _reference_time_point = -1;
    }

    /// Update
    void update();

    /// Set floating time point
    void set_floating_time_point(int floating_time_point) { _floating_time_point = floating_time_point; }

    /// Set reference time point
    void set_reference_time_point(int reference_time_point) { _reference_time_point = reference_time_point; }

protected:

    /// Check parameters
    virtual void check_parameters();

    /// Parse parameter file
    virtual void parse_parameter_file();

    /// Registration object
    std::shared_ptr<reg_f3d<T> > _registration_sptr;

    /// Floating time point
    int _floating_time_point;
    /// Reference time point
    int _reference_time_point;
};

#endif
