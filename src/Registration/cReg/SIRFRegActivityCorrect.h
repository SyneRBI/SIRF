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
\brief Resampling class based on nifty resample

\author Richard Brown
\author CCP PETMR
*/

#ifndef _SIRFREGACTIVITYCORRECT_H_
#define _SIRFREGACTIVITYCORRECT_H_

#include <nifti1_io.h>
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

/// Activity correct an image
class SIRFRegActivityCorrect
{
public:

    /// Constructor
    SIRFRegActivityCorrect();

    /// Destructor
    virtual ~SIRFRegActivityCorrect() {}

    /// Set input .nii image
    void set_input(const nifti_image *input_image_ptr)
    {
        _input_image_sptr     = std::make_shared<nifti_image>(*input_image_ptr);
        _input_image_filename = "";
    }

    /// Set input .nii filename
    void set_input_image_filename(const std::string input_image_filename)
    {
        _input_image_filename = input_image_filename;
        _input_image_sptr.reset();
    }

    /// Set initial activity (Bq)
    void set_initial_activity(float initial_activity) { _initial_activity = initial_activity; }

    /// Set half life (s)
    void set_half_life(float half_life) { _half_life = half_life; }

    /// Set start (s) (time after initial activity that image starts)
    void set_start(float start) { _start = start; }

    /// Set stop (s) (time after initial activity that image stops)
    void set_stop(float stop) { _stop = stop; }

    /// Update
    void update();

    /// Get output
    nifti_image *get_output() const { return _output_image_sptr.get(); }

    /// Save output image to file
    void save_output(const std::string filename) const;

protected:

    /// Check parameters
    virtual void check_parameters();

    /// Input image filename
    boost::filesystem::path      _input_image_filename;
    /// Input image
    std::shared_ptr<nifti_image> _input_image_sptr;

    /// Initial activitiy (s)
    float _initial_activity;
    /// Half-life (s)
    float _half_life;
    /// Start time (s)
    float _start;
    /// Stop time (s)
    float _stop;

    /// Output image
    std::shared_ptr<nifti_image> _output_image_sptr;
};

#endif
