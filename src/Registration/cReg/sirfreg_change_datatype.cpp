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
\brief Convert affine transformation matrix to displacement or deformation field image(s)

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include <vector>
#include <NiftiImage.h>

using namespace std;
using namespace sirf;

/// Print usage
void print_usage(const vector<string> &datatypes)
{
    cout << "\n\n\n*** Usage: sirfreg_change_datatype output_filename input_filename desired_datatype ***\n\n";
    cout << "Supported datatypes:\n";
    for (unsigned i=0; i<datatypes.size(); ++i)
        cout << "\t" << datatypes[i] << "\n";
    cout << "\n\n";
}

/// main
int main(int argc, char* argv[])
{

    try {
        // Supported datatypes
        vector<string> datatypes;
        datatypes.push_back("signed short");
        datatypes.push_back("signed int");
        datatypes.push_back("float");
        datatypes.push_back("double");
        datatypes.push_back("unsigned char");
        datatypes.push_back("unsigned short");
        datatypes.push_back("unsigned int");
        datatypes.push_back("signed long long");
        datatypes.push_back("unsigned long long");
        datatypes.push_back("long double");

        if (argc != 3) {
            print_usage(datatypes);
            return EXIT_SUCCESS;
        }

        // Read image
        NiftiImage im(argv[2]);

        // Change datatype
        if      (strcmp(argv[3], "signed short")       == 0) im.change_datatype<signed short>();
        else if (strcmp(argv[3], "signed int")         == 0) im.change_datatype<signed int>();
        else if (strcmp(argv[3], "float")              == 0) im.change_datatype<float>();
        else if (strcmp(argv[3], "double")             == 0) im.change_datatype<double>();
        else if (strcmp(argv[3], "unsigned char")      == 0) im.change_datatype<unsigned char>();
        else if (strcmp(argv[3], "unsigned short")     == 0) im.change_datatype<unsigned short>();
        else if (strcmp(argv[3], "unsigned int")       == 0) im.change_datatype<unsigned int>();
        else if (strcmp(argv[3], "signed long long")   == 0) im.change_datatype<signed long long>();
        else if (strcmp(argv[3], "unsigned long long") == 0) im.change_datatype<unsigned long long>();
        else if (strcmp(argv[3], "long double")        == 0) im.change_datatype<long double>();
        else {
            print_usage(datatypes);
            return EXIT_SUCCESS;
        }

        // Save output
        im.save_to_file(argv[1]);

    // If there was an error
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
