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
\brief Perform weighted mean of multiple images.

\author Richard Brown
\author CCP PETMR
*/

#pragma once

#include <vector>
#include <memory>

namespace sirf {

// Forward delcarations
template<class dataType> class NiftiImageData;

/// Calculate the weighted mean of a set of images
template<class dataType>
class ImageWeightedMean
{
public:

    /// Constructor
    ImageWeightedMean();

    /// Destructor
    ~ImageWeightedMean() {}

    /// Add an image (from NiftImage) and its corresponding weight
    void add_image(const NiftiImageData<dataType> &image, const float weight);

    /// Process
    void process();

    /// Get output
    const std::shared_ptr<const NiftiImageData<dataType> > get_output_sptr() const { return _output_image_sptr; }

protected:

    /// Check if its possible to calculate the mean
    void check_can_do_mean() const;

    /// Bool to check if update is necessary
    bool _need_to_update;
    /// Vector of input images
    std::vector<std::shared_ptr<const NiftiImageData<dataType> > > _input_image_sptrs;
    /// Vector of weights
    std::vector<float> _weights;
    /// Output image
    std::shared_ptr<NiftiImageData<dataType> > _output_image_sptr;

};
}
