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

SIRFRegImageWeightedMean::SIRFRegImageWeightedMean()
{
    _need_to_update = true;
}

void SIRFRegImageWeightedMean::add_image(const string &filename, const float weight)
{
    // Open image
    std::shared_ptr<nifti_image> input_image_sptr;
    SIRFRegMisc::open_nifti_image(input_image_sptr, filename);

    // Use other function to add image to list of vectors
    this->add_image(input_image_sptr.get(), weight);
}

void SIRFRegImageWeightedMean::add_image(const SIRFImageData &image, const float weight)
{
    // Add image to vector
    _input_images.push_back(image);
    _weights.push_back(weight);

    _need_to_update = true;
}

void SIRFRegImageWeightedMean::update()
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
    _output_image = _input_images[0];

    // Get the data and cast it to float
    float *output_data_ptr = static_cast<float*>(_output_image.get_image_as_nifti()->data);

    // Set all of the output image's voxels to 0
    for (unsigned i=0; i<_output_image.get_image_as_nifti()->nvox; i++) output_data_ptr[i] = 0.;

    // Loop over each input image
    for (unsigned i=0; i<_input_images.size(); i++) {

        // Get the data and cast it to float for the ith input image
        float *input_data_ptr = static_cast<float*>(_input_images[i].get_image_as_nifti()->data);

        // Loop over each voxel
        for (unsigned j=0; j<_output_image.get_image_as_nifti()->nvox; j++) {

            // Add in the weighted contribution of the jth voxel of the ith image
            output_data_ptr[j] += input_data_ptr[j] * normalised_weights[i];
        }
    }

    // Once the update is done, set the need_to_update flag to false
    _need_to_update = false;
}

void SIRFRegImageWeightedMean::save_image_to_file(const string &filename) const
{
    _output_image.save_to_file(filename);
}

void SIRFRegImageWeightedMean::check_can_do_mean() const
{
    std::string errors;

    // Check that num_images > 0. If not, throw error
    if (_input_images.size() == 0)
        errors += "Need to add images to be able to do weighted mean.\n";

    // Print all the info
    std::cout << "\n\nChecking that images can be combined...\n\n";
    std::cout << "im || datatype | ndim | nvox    | nx  | ny  | nz | dx     | dy     | dz  \n";
    std::cout << "---++----------+------+---------+-----+-----+----+--------+--------+-----\n";
    for (unsigned i=0; i<_input_images.size(); i++) {
        std::cout << i << "  || ";
        std::string data_type = nifti_datatype_string(_input_images[i].get_image_as_nifti()->datatype);
        std::cout << data_type             <<   "  | ";
        std::cout << _input_images[i].get_image_as_nifti()->ndim << "    | ";
        std::cout << _input_images[i].get_image_as_nifti()->nvox <<    " | ";
        std::cout << _input_images[i].get_image_as_nifti()->nx   <<    " | ";
        std::cout << _input_images[i].get_image_as_nifti()->ny   <<    " | ";
        std::cout << _input_images[i].get_image_as_nifti()->nz   <<    " | ";
        std::cout << _input_images[i].get_image_as_nifti()->dx   <<    " | ";
        std::cout << _input_images[i].get_image_as_nifti()->dy   <<    " | ";
        std::cout << _input_images[i].get_image_as_nifti()->dz   << "\n";
    }
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
