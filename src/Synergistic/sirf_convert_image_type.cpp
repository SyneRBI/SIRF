/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2019 University College London

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
\ingroup Synergistic
\brief Convert any SIRF image type into any other.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"

using namespace sirf;

static const std::shared_ptr<ImageData> image_as_sptr(const std::string &filename, const std::string &engine, const bool verbose)
{
    std::shared_ptr<ImageData> img_sptr;

    if (strcmp(engine.c_str(), "Reg") == 0) {
        std::shared_ptr<NiftiImageData<float> > nifti_sptr = std::make_shared<NiftiImageData3D<float> >(filename);
        if (verbose) nifti_sptr->print_header();
        img_sptr = nifti_sptr;
    }
    else if (strcmp(engine.c_str(), "STIR") == 0) {
        img_sptr = std::make_shared<STIRImageData>(filename);
    }

    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        std::shared_ptr<GadgetronImagesVector> gadgetron_sptr(new GadgetronImagesVector);
		gadgetron_sptr->read(filename);
        if (verbose) gadgetron_sptr->print_header(0);
        img_sptr = gadgetron_sptr;
    }
    else
        throw std::runtime_error("unknown engine - " + engine + ".\n");

    // If verbose print geom info
    if (verbose) img_sptr->get_geom_info_sptr()->print_info();

    // return
    return img_sptr;
}

static void convert_and_write_image(const std::string &filename, const std::string &engine, const std::shared_ptr<ImageData> &in_img_sptr, const std::string &param_file, const bool verbose)
{
    if (strcmp(engine.c_str(), "Reg") == 0) {
        NiftiImageData<float> im(*in_img_sptr);
        im.write(filename);
        if (verbose) {
            im.print_header();
            im.get_geom_info_sptr()->print_info();
        }
    }
    else if (strcmp(engine.c_str(), "STIR") == 0) {
        STIRImageData im(*in_img_sptr);
        if (param_file.empty())
            im.write(filename);
        else
            im.write(filename,param_file);
        if (verbose)
            im.get_geom_info_sptr()->print_info();
    }
    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        throw std::runtime_error("Converter to GadgetronImagesVector not yet implemented.\n");
    }
    else
        throw std::runtime_error("unknown engine - " + engine + ".\n");
}

/// print usage
void print_usage(const std::string &app_name)
{
    std::cout << "\n*** " << app_name << " [-h/--help][-v/--verbose] out_filename out_engine in_filename in_engine [STIR_param_file]***\n";

    // Required arguments
    std::cout << "\n  Required arguments:\n";
    std::cout << "    out_filename:\t\toutput image filename (with or without extension)\n";
    std::cout << "    out_engine:\t\t\tengine for output image (Reg/Nifti/nii, STIR or Gadgetron/ISMRMRD/h5)\n";
    std::cout << "    in_filename:\t\tinput image filename\n";
    std::cout << "    in_engine:\t\t\tengine for reading input image (Reg/Nifti/nii, STIR or Gadgetron/ISMRMRD/h5)\n";

    // Optional arguments
    std::cout << "\n  Optional arguments:\n";
    std::cout << "    STIR_param_file:\t\tif the out_engine is set to STIR you can optionally "
                 "supply an output parameter file (such as /examples/parameter_files/STIR_output_file_format_nifti.par)\n";

    // Optional flags
    std::cout << "\n  Optional flags:\n";
    std::cout << "    -h/--help:\t\t\tprint this help message\n";
    std::cout << "    -v/--verbose:\t\t\tprint headers of input and output images where possible\n";
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

        // Check for help
        for (unsigned i=1; i<unsigned(argc); ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
            }
        }

        // Check for verbose
        bool verbose = false;
        if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--verbose") == 0) {
            verbose = true;
            ++argv;
            --argc;
        }

        // Check number of inputs
        if(argc < 5 && argc < 7) {
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
        }

        // Get filenames and engines
        std::string out_filename = argv[1];
        std::string out_engine   = argv[2];
        std::string in_filename  = argv[3];
        std::string in_engine    = argv[4];

        // For engines, convert Nifti and nii to Reg and convert ISMRMRD and h5 to Gadgetron
        if (strcmp(in_engine.c_str(), "Nifti") == 0 || strcmp(in_engine.c_str(), "nii") == 0)
            in_engine = "Reg";
        if (strcmp(out_engine.c_str(), "Nifti") == 0 || strcmp(out_engine.c_str(), "nii") == 0)
            out_engine = "Reg";
        if (strcmp(in_engine.c_str(), "ISMRMRD") == 0 || strcmp(in_engine.c_str(), "h5") == 0)
            in_engine = "Gadgetron";
        if (strcmp(out_engine.c_str(), "ISMRMRD") == 0 || strcmp(out_engine.c_str(), "h5") == 0)
            out_engine = "Gadgetron";

        // If there's a parameter file, check that the output engine is STIR.
        // Else throw error
        std::string param_file = "";
        if (argc == 6) {
            if (strcmp(out_engine.c_str(), "STIR") == 0)
                param_file = argv[5];
            else {
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
        }

        // Read input image
        std::shared_ptr<ImageData> in_img_sptr = image_as_sptr(in_filename, in_engine, verbose);

        // Convert and write
        convert_and_write_image(out_filename, out_engine, in_img_sptr, param_file, verbose);
    }

    // If there was an error
    catch(const std::exception &error) {
        std::cerr << "\nError encountered:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
