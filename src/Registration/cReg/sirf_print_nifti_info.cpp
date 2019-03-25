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
\brief Dump info of one or more nifti images

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData.h"
#include "vector"

using namespace sirf;

/// main
int main(int argc, char* argv[])
{
    try {
        if (argc < 2) {
            std::cout << "\nsirf_print_nifti_info filename1 [filename2 [filename3 [...]]]\n";
            return EXIT_SUCCESS;
        }

        // Vector of images
        int num_images = argc - 1;
        std::vector<std::shared_ptr<const NiftiImageData<float> > > ims;
        std::vector<const NiftiImageData<float>*> ims_ptr;

        // Read all the images
        for (int i=0; i<num_images; ++i) {
            ims.push_back(std::make_shared<const NiftiImageData<float> >(argv[i+1]));
            ims_ptr.push_back(ims[i].get());
        }

        // Print info
        NiftiImageData<float>::print_headers(ims_ptr);

    // If there was an error
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
