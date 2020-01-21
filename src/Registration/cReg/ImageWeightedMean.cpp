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

#include "sirf/Reg/ImageWeightedMean.h"
#include "sirf/Reg/NiftiImageData.h"
#include "sirf/Reg/NiftyResample.h"
#include <iostream>

using namespace sirf;

template<class dataType>
ImageWeightedMean<dataType>::ImageWeightedMean()
{
    _need_to_update = true;
}

template<class dataType>
void ImageWeightedMean<dataType>::add_image(const NiftiImageData<dataType> &image, const float weight)
{
    // Add image to vector
    _input_image_sptrs.push_back(std::make_shared<const NiftiImageData<dataType> >(image));
    _weights.push_back(weight);

    _need_to_update = true;
}

template<class dataType>
void ImageWeightedMean<dataType>::process()
{
    // Only process if you need to
    if (!_need_to_update) return;

    // Do all the various checks before performing average
    check_can_do_mean();

    // Check the image metdata. If they don't match, resample
    std::vector<std::shared_ptr<const NiftiImageData<dataType> > > images_sptr;
    for (unsigned i=0; i<_input_image_sptrs.size(); i++) {
        if (NiftiImageData<dataType>::do_nifti_image_metadata_match(*_input_image_sptrs[0],*_input_image_sptrs[i], false))
            images_sptr.push_back(_input_image_sptrs[i]);
        else {
            std::cout << "\nImageWeightedMean: Mismatch in metadata between images 0 and " << i << ", so resampling taking first image as reference...\n";
            NiftyResample<float> resample;
            resample.set_interpolation_type_to_nearest_neighbour();
            resample.set_reference_image(_input_image_sptrs[0]);
            resample.set_floating_image(_input_image_sptrs[i]);
            images_sptr.push_back(
                        std::dynamic_pointer_cast<const NiftiImageData<dataType> >(
                            resample.forward(_input_image_sptrs[i])));
        }
    }

    // Need to normalise the weights so that sum = 1
    float sum_of_weights = 0.;
    std::vector<float> normalised_weights = _weights;
    for (unsigned i=0; i<_weights.size(); i++) sum_of_weights += _weights[i];
    for (unsigned i=0; i<_weights.size(); i++) normalised_weights[i] /= sum_of_weights;

    // Create a copy of the first image to use as a template for the output
    _output_image_sptr = std::make_shared<NiftiImageData<dataType> >(*_input_image_sptrs[0]);

    // Set all of the output image's voxels to 0
    _output_image_sptr->fill(0.F);

    // Loop over each input image and each voxel
    for (unsigned i=0; i<_input_image_sptrs.size(); i++)
        for (int j=0; j<int(_output_image_sptr->get_raw_nifti_sptr()->nvox); j++)
            // Add in the weighted contribution of the jth voxel of the ith image
            (*_output_image_sptr)(j) += (*_input_image_sptrs[i])(j) * normalised_weights[i];

    // Once the processing is done, set the need_to_update flag to false
    _need_to_update = false;
}

template<class dataType>
void ImageWeightedMean<dataType>::check_can_do_mean() const
{
    // Check that num_images > 0. If not, throw error
    if (_input_image_sptrs.size() == 0)
        throw std::runtime_error("Need to add images to be able to do weighted mean.");

    std::cout << "\nAll images match, we can calculate their weighted average.\n";
}

namespace sirf {
template class ImageWeightedMean<float>;
}
