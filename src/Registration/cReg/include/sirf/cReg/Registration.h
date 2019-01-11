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
\brief Base class for all SIRF registration.

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include <string>
#include <vector>
#include <memory>

namespace sirf {

/// Forward declarations
template<class dataType> class Transformation;
class ImageData;

/*!
\ingroup Registration
\brief Base class for all SIRF registration.

Registration classes can be set up with parameter files or inputting arguments directly with the set_parameter() method.

If using a parameter file, it should have interfile-like syntax.
The variables will be stored as a vector of floats and converted into the required type (int, unsigned int, etc) if necessary.
Multiple variables for a given parameter should be comma separated.
Spaces and tabs will be ignored.
For the title, it doesn't matter what is written as it will be ignored, but something has to be there (otherwise the first parameter will be ignored).
Possible parameters are all the Set<something> methods for each class (e.g., nifty_aladin::SetPerformRigid) and should be written in the parameter file without the "Set" (e.g., PerformRigid).

An example is given below:
\verbatim
    SomeTitle :=
        ReferenceTimePoint := 1
        FloatingTimePoint := 2
        LinearEnergyWeights := 1.5,1
        AdditiveMC :=
    end :=
\endverbatim

More examples can be found in data/examples/Registration/paramFiles

\author Richard Brown
\author CCP PETMR
*/
template<class dataType>
class Registration
{
public:

    /// Constructor
    Registration() {}

    /// Destructor
    virtual ~Registration() {}

    /// Set parameter file
    void set_parameter_file(const std::string &parameter_filename) { _parameter_filename = parameter_filename; }

    /// Set reference image
    void set_reference_image(const std::shared_ptr<const ImageData> reference_image_sptr) { _reference_image_sptr = reference_image_sptr; }

    /// Set floating image
    void set_floating_image(const std::shared_ptr<const ImageData> floating_image_sptr) { _floating_image_sptr = floating_image_sptr; }

    /// Process
    virtual void process() = 0;

    /// Get registered image
    const std::shared_ptr<const ImageData> get_output_sptr() const { return _warped_image_sptr; }

    /// Get forward deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_forward_sptr() const = 0;

    /// Get inverse deformation field image
    virtual const std::shared_ptr<const Transformation<dataType> > get_deformation_field_inverse_sptr() const = 0;

    /// Get forward displacement field image
    const std::shared_ptr<const Transformation<dataType> > get_displacement_field_forward_sptr() const { return _disp_image_forward_sptr; }

    /// Get inverse displacement field image
    const std::shared_ptr<const Transformation<dataType> > get_displacement_field_inverse_sptr() const { return _disp_image_inverse_sptr; }

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

    /// Check parameters
    virtual void check_parameters() const;

    /// Set any extra parameters
    virtual void set_parameters() = 0;

    /// Store extra parameters. Only apply them after parsing.
    std::vector<std::string> _extra_params;

    /// Parameter filename
    std::string _parameter_filename;

    /// Reference image
    std::shared_ptr<const ImageData> _reference_image_sptr;
    /// Floating image
    std::shared_ptr<const ImageData> _floating_image_sptr;
    /// Warped image
    std::shared_ptr<ImageData> _warped_image_sptr;

    /// Forward displacement field image
    std::shared_ptr<Transformation<dataType> > _disp_image_forward_sptr;
    /// Inverse displacement field image
    std::shared_ptr<Transformation<dataType> > _disp_image_inverse_sptr;

    /// Floating mask
    std::shared_ptr<const ImageData> _floating_mask_sptr;
    /// Reference mask
    std::shared_ptr<const ImageData> _reference_mask_sptr;
};
}
