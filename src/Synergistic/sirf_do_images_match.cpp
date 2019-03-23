/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2019 Rutherford Appleton Laboratory STFC

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
\brief Check if two SIRF images match to a given precision.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"


using namespace sirf;

static std::shared_ptr<const NiftiImageData3D<float> > image_as_sptr(const std::string &filename, const std::string &engine)
{
    if (strcmp(engine.c_str(), "Nifti") == 0)
        return std::make_shared<const NiftiImageData3D<float> >(filename);
    else if (strcmp(engine.c_str(), "STIR") == 0)
        return std::make_shared<const NiftiImageData3D<float> >(STIRImageData(filename));
    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        std::shared_ptr<GadgetronImageData> sptr_img(new GadgetronImagesVector);
		sptr_img->read(filename);
        return std::make_shared<const NiftiImageData3D<float> >(*sptr_img);
    }
    else
        throw std::runtime_error("unknown engine - " + engine + ".\n");
}

/// print usage
void print_usage()
{
    std::cout << "\n*** sirf_do_images_match usage ***\n";

    // Required flags
    std::cout << "\n  Required flags:\n";
    std::cout << "    -im1:\t\timage 1\n";
    std::cout << "    -im2:\t\timage 2\n";

    // Optional flags
    std::cout << "\n  Optional flags:\n";
    std::cout << "    -eng_im1:\t\tengine to open image1\n";
    std::cout << "    -eng_im2:\t\tengine to open image 2\n";
    std::cout << "    -accu:\t\tAccuracy (compared to avg. of max of two images)\n";
}

/// throw error
void err(const std::string message)
{
    std::cerr << "\n" << message << "\n";
    exit(EXIT_FAILURE);
}

/// main
int main(int argc, char* argv[])
{
    try {

        std::string im1_filename = "", im2_filename = "";
        std::string eng_im1 = "Nifti", eng_im2 = "Nifti";
        float accuracy = 1.E-3F;

        // Loop over all input arguments (ignore first argument (name of executable))
        argc--; argv++;
        while (argc>0) {

            // help
            if (strcmp(argv[0], "-h") == 0) {
                print_usage();
                exit(EXIT_SUCCESS);
            }

            // im1
            if (strcmp(argv[0], "-im1") == 0) {
                if (argc<2)
                    err("Option '-im1' expects a (filename) argument.");
                im1_filename = argv[1];
                argc-=2; argv+=2;
            }
            // im2
            else if (strcmp(argv[0], "-im2") == 0) {
                if (argc<2)
                    err("Option '-im2' expects a (filename) argument.");
                im2_filename = argv[1];
                argc-=2; argv+=2;
            }
            // im1 engine
            else if (strcmp(argv[0], "-eng_im1") == 0) {
                if (argc<2)
                    err("Option '-eng_im1' expects a (string) argument.");
                eng_im1 = argv[1];
                argc-=2; argv+=2;
            }
            // im2 engine
            else if (strcmp(argv[0], "-eng_im2") == 0) {
                if (argc<2)
                    err("Option '-eng_im2' expects a (string) argument.");
                eng_im2 = argv[1];
                argc-=2; argv+=2;
            }
            // accuracy
            else if (strcmp(argv[0], "-accu") == 0) {
                if (argc<2)
                    err("Option '-accu' expects a (numerical) argument.");
                accuracy = float(atof(argv[1]));
                argc-=2; argv+=2;
            }

            // Unknown argument
            else
                err("Unknown option '" + std::string(argv[0]) + "'.");
        }

        // Check filenames aren't blank
        if (im1_filename.size() == 0)
            err("Error: -im1 required.");
        if (im2_filename.size() == 0)
            err("Error: -im2 required.");

        // Get images as NiftiImages
        std::shared_ptr<const NiftiImageData3D<float> > im1 = image_as_sptr(im1_filename,eng_im1);
        std::shared_ptr<const NiftiImageData3D<float> > im2 = image_as_sptr(im2_filename,eng_im2);

        // Compare
        bool result = NiftiImageData3D<float>::are_equal_to_given_accuracy(*im1,*im2,accuracy);

        if (result)
            std::cout << "\nImages match!\n";
        else
            std::cout << "\nImages don't match.\n";
    }

    // If there was an error
    catch(const std::exception &error) {
        std::cerr << "\nError encountered:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
