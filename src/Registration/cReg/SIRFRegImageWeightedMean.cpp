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
\brief Perform weighted mean of multiple images.

\author Richard Brown
\author CCP PETMR
*/

#include "SIRFRegImageWeightedMean.h"
#include "SIRFRegMisc.h"
#include <boost/filesystem.hpp>
#include <iostream>

using namespace std;
using namespace sirf;

template<class ImType>
SIRFRegImageWeightedMean<ImType>::SIRFRegImageWeightedMean()
{
    _need_to_update = true;
}

template<class ImType>
void SIRFRegImageWeightedMean<ImType>::add_image(const std::string &filename, const float weight)
{
    // Use other function to add image to list of vectors
    this->add_image(ImType(filename), weight);
}

template<class ImType>
void SIRFRegImageWeightedMean<ImType>::add_image(const ImType &image, const float weight)
{
    // Add image to vector
    _input_images.push_back(image.deep_copy());
    _weights.push_back(weight);

    _need_to_update = true;
}

template<class ImType>
void SIRFRegImageWeightedMean<ImType>::update()
{
    // Only update if you need to
    if (!_need_to_update) return;

    // Do all the various checks before performing average
    check_can_do_mean();

    // Need to normalise the weights so that sum = 1
    float sum_of_weights = 0.;
    vector<float> normalised_weights = _weights;
    for (unsigned i=0; i<_weights.size(); i++) sum_of_weights += _weights[i];
    for (unsigned i=0; i<_weights.size(); i++) normalised_weights[i] /= sum_of_weights;

    // Create a copy of the first image to use as a template for the output
    _output_image = _input_images[0].deep_copy();

    // Change to double to minimise rounding errors. Get the data.
    SIRFRegMisc::change_datatype<double>(_output_image);
    double *output_data_ptr = static_cast<double*>(_output_image.get_image_as_nifti()->data);

    // Set all of the output image's voxels to 0
    _output_image.fill(0.F);

    // Loop over each input image
    for (unsigned i=0; i<_input_images.size(); i++) {

        // Create a temporary copy of the image so that we can change the datatype
        ImType temp = _input_images[i].deep_copy();
        SIRFRegMisc::change_datatype<double>(temp);

        // Get the data and cast it to float for the ith input image
        double *input_data_ptr = static_cast<double*>(temp.get_image_as_nifti()->data);

        // Loop over each voxel
        for (unsigned j=0; j<_output_image.get_image_as_nifti()->nvox; j++) {

            // Add in the weighted contribution of the jth voxel of the ith image
            output_data_ptr[j] += input_data_ptr[j] * double(normalised_weights[i]);
        }
    }

    // Put the output type back so that it matches the input type
    if (_input_images[0].get_image_as_nifti()->datatype == DT_INT16)   SIRFRegMisc::change_datatype<signed short>  (_output_image);
    if (_input_images[0].get_image_as_nifti()->datatype == DT_INT32)   SIRFRegMisc::change_datatype<signed int>    (_output_image);
    if (_input_images[0].get_image_as_nifti()->datatype == DT_FLOAT32) SIRFRegMisc::change_datatype<float>         (_output_image);
    if (_input_images[0].get_image_as_nifti()->datatype == DT_FLOAT64) SIRFRegMisc::change_datatype<double>        (_output_image);
    if (_input_images[0].get_image_as_nifti()->datatype == DT_UINT8)   SIRFRegMisc::change_datatype<unsigned char> (_output_image);
    if (_input_images[0].get_image_as_nifti()->datatype == DT_UINT16)  SIRFRegMisc::change_datatype<unsigned short>(_output_image);
    if (_input_images[0].get_image_as_nifti()->datatype == DT_UINT32)  SIRFRegMisc::change_datatype<unsigned int>  (_output_image);

    // Once the update is done, set the need_to_update flag to false
    _need_to_update = false;
}

template<class ImType>
void SIRFRegImageWeightedMean<ImType>::save_image_to_file(const string &filename) const
{
    _output_image.save_to_file(filename);
}

template<class ImType>
void SIRFRegImageWeightedMean<ImType>::check_can_do_mean() const
{
    // Check that num_images > 0. If not, throw error
    if (_input_images.size() == 0)
        throw std::runtime_error("Need to add images to be able to do weighted mean.");

    bool can_do_mean = true;

    // Check each of the images against all the other images
    for (unsigned i=0; i<_input_images.size(); i++) {
        for (unsigned j=i+1; j<_input_images.size(); j++) {

            std::cout << "\nComparing input images " << i << " and " << j << "...\n";
            if (!SIRFRegMisc::do_nifti_image_metadata_match(_input_images[i],_input_images[j])) can_do_mean = false;
        }
    }

    // if there were any errors, throw them
    if (!can_do_mean) {
        throw std::runtime_error("There is a mismatch in images. Cannot calculate their weighted mean.");
    }

    std::cout << "\nAll images match, we can calculate their weighted average.\n";
}

namespace sirf {
// Put the instantiations of the template class at the END of the file!
template class SIRFRegImageWeightedMean<SIRFImageData>;
template class SIRFRegImageWeightedMean<SIRFImageDataDeformation>;
}
