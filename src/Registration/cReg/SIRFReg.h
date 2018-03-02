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
\brief Base class for all SIRF registration.

The parameter file should have interfile-like syntax.
The variables will be stored as a vector of floats and converted into the required type (int, unsigned int, etc) if necessary.
Multiple variables for a given parameter should be comma separated.
Spaces and tabs will be ignored.
For the title, it doesn't matter what is written as it will be ignored, but something has to be there (otherwise the first parameter will be ignored).
Possible parameters are all the Set<something> methods for each class (e.g., nifty_aladin::SetPerformRigid) and should be written in the parameter file without the "Set" (e.g., PerformRigid).

An example is given below:
    SomeTitle :=
        ReferenceTimePoint := 1
        FloatingTimePoint := 2
        LinearEnergyWeights := 1.5,1
        AdditiveMC :=
    end :=

More examples can be found in // Need to give path

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFREG_H_
#define _SIRFREG_H_

#include <stdexcept>
#include <nifti1_io.h>
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

class SIRFReg
{
public:

    /// Constructor
    SIRFReg() {}

    /// Destructor
    virtual ~SIRFReg() {}

    /// Set parameter file
    void set_parameter_file(const std::string parameter_filename) { _parameter_filename = parameter_filename; }

    /// Set reference image filename
    void set_reference_image_filename(const std::string filename)
    {
        _reference_image_filename = filename;
        _reference_image_sptr.reset();
    }

    /// Set reference image
    void set_reference_image(const nifti_image *reference_image)
    {
        _reference_image_sptr     = std::make_shared<nifti_image>(*reference_image);
        _reference_image_filename = "";
    }

    /// Set floating image filename
    void set_floating_image_filename(const std::string filename)
    {
        _floating_image_filename = filename;
        _floating_image_sptr.reset();
    }

    /// Set floating image
    void set_floating_image(const nifti_image *floating_image)
    {
        _floating_image_sptr = std::make_shared<nifti_image>(*floating_image);
        _floating_image_filename = "";
    }

    /// Update
    virtual void update() = 0;

    /// Get registered image
    nifti_image *get_output() const { return _warped_image_sptr.get(); }

    /// Save warped image to file
    void save_warped_image(const std::string filename) const;

protected:

    /// Parse parameter file
    virtual void parse_parameter_file() = 0;

    /// Check parameters
    virtual void check_parameters();

    boost::filesystem::path      _parameter_filename;
    boost::filesystem::path      _floating_image_filename;
    boost::filesystem::path      _reference_image_filename;
    std::shared_ptr<nifti_image> _floating_image_sptr;
    std::shared_ptr<nifti_image> _reference_image_sptr;
    std::shared_ptr<nifti_image> _warped_image_sptr;
};

#endif
