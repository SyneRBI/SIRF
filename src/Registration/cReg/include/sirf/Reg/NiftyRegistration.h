/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2017 - 2020 University College London

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
\brief Base class for all NiftyReg registrations.

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include "sirf/Reg/NiftiBasedRegistration.h"

namespace sirf {

/// Forward declarations
template<class dataType> class NiftiImageData3D;

/*!
\ingroup Registration
\brief Base class for all NiftyReg registrations.

\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class NiftyRegistration : public NiftiBasedRegistration<dataType>
{
public:

    /// Constructor
    NiftyRegistration();

    /// Destructor
    virtual ~NiftyRegistration() {}

    /// Set parameter file
    void set_parameter_file(const std::string &parameter_filename) { _parameter_filename = parameter_filename; }

    /// Set string parameter. Check if any set methods match the method given by par.
    /// If so, set the value given by arg. Convert to float/int etc., as necessary.
    /// Up to 2 arguments, leave blank if unneeded. These are applied after parsing
    /// the parameter file.
    void set_parameter(const std::string &par, const std::string &arg1 = "", const std::string &arg2 = "");

    /// Set reference mask
    void set_reference_mask(const std::shared_ptr<const ImageData> reference_mask_sptr) { _reference_mask_sptr = reference_mask_sptr; }

    /// Set floating mask
    void set_floating_mask(const std::shared_ptr<const ImageData> floating_mask_sptr)   {  _floating_mask_sptr = floating_mask_sptr;  }

protected:

    /// Parse parameter file
    virtual void parse_parameter_file() = 0;

    /// Set up inputs
    void set_up_inputs();

    /// Set any extra parameters
    virtual void set_parameters() = 0;

    /// Store extra parameters. Only apply them after parsing.
    std::vector<std::string> _extra_params;

    /// Parameter filename
    std::string _parameter_filename;

    /// Floating mask
    std::shared_ptr<const ImageData> _floating_mask_sptr;
    /// Reference mask
    std::shared_ptr<const ImageData> _reference_mask_sptr;

    /// Floating mask (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _floating_mask_nifti_sptr;
    /// Reference mask (as NiftiImageData3D)
    std::shared_ptr<const NiftiImageData3D<dataType> > _reference_mask_nifti_sptr;
};
}
