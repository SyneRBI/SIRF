/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccpsynerbi.ac.uk/).

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
\ingroup Gadgetron Extensions
\brief Auxiliary functions for MR related C++ tests.

\author Johannes Mayer
*/

#include <iostream>
#include "mrtest_auxiliary_funs.h"

int main ( int argc, char* argv[])
{
    try{    

        if(argc != 3)
            throw std::runtime_error("Please supply exactly two arguments: full path to input and output files.");
        
        const char * path_in  = argv[1];
        const char * path_out = argv[2];
                
        sirf::set_acq_default_orientation(path_in, path_out);

        return 0;
    }
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
