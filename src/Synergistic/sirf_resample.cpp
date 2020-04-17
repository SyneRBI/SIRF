/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2019 University College London

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
\brief Perform resampling with any type of SIRF image and any SIRF resampling algorithm.

If no transformations are given, identity will be used.
If multiple transformations are given, they will be applied in the order they were given.

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Reg/NiftiImageData3DDeformation.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include "sirf/STIR/stir_data_containers.h"


using namespace sirf;

static std::shared_ptr<const ImageData> image_as_sptr(const std::string &filename, const std::string &engine)
{
    if (strcmp(engine.c_str(), "Nifti") == 0)
        return std::make_shared<const NiftiImageData<float> >(filename);
    else if (strcmp(engine.c_str(), "STIR") == 0)
        return std::make_shared<const STIRImageData>(filename);
    else if (strcmp(engine.c_str(), "Gadgetron") == 0) {
        std::shared_ptr<GadgetronImageData> sptr_img(new GadgetronImagesVector);
		sptr_img->read(filename);
        return sptr_img;
    }
    else
        throw std::runtime_error("unknown engine - " + engine + ".\n");
}

static std::shared_ptr<Resample<float> > algo_as_sptr(const std::string &algorithm)
{
    std::cout << "\nUsing " << algorithm << " resampling algorithm...\n";
    if      (strcmp(algorithm.c_str(), "niftyreg") == 0)
        return std::make_shared<NiftyResample<float> >();
    else
        throw std::runtime_error("Synergistic_registration: unknown algorithm - " + algorithm + ".\n");
}

/// print usage
void print_usage()
{
    std::cout << "\n*** sirf_resample usage ***\n";

    // Required flags
    std::cout << "\n  Required flags:\n";
    std::cout << "    -algo:\t\tResampling algorithm (currently only niftyreg)\n";
    std::cout << "    -ref:\t\treference image\n";
    std::cout << "    -flo:\t\tfloating image\n";

    // Optional flags
    std::cout << "\n  Optional flags:\n";
    std::cout << "    -eng_ref:\t\tengine to open reference image (Nifti, STIR, Gadgetron)\n";
    std::cout << "    -eng_flo:\t\tengine to open floating image (Nifti, STIR, Gadgetron)\n";
    std::cout << "    -output:\t\toutput image filename\n";
    std::cout << "    -interp:\t\tinterpolation (0=NN, 1=linear, 3=cubic, 4=spline)\n";
    std::cout << "    -add_affine:\tadd affine transformation\n";
    std::cout << "    -add_def:\t\tadd deformation transformation\n";
    std::cout << "    -add_disp:\t\tadd displacement transformation\n";
    std::cout << "    -adj:\t\tadjoint transformation. Give ref and flo as you would in the forward case.\n";

}

/// throw error
[[ noreturn ]] void err(const std::string message)
{
    std::cerr << "\n" << message << "\n";
    exit(EXIT_FAILURE);
}

/// main
int main(int argc, char* argv[])
{
    try {

        std::string ref_filename = "", flo_filename = "";
        std::string eng_ref = "Nifti", eng_flo = "Nifti";
        std::string output = "output";
        std::vector<std::shared_ptr<const Transformation<float> > > trans;
        std::string algo = "niftyreg";
        Resample<float>::InterpolationType interp = Resample<float>::NEARESTNEIGHBOUR;
        float pad = 0;
        bool pad_set = false;
        bool forward = true;

        // Loop over all input arguments (ignore first argument (name of executable))
        argc--; argv++;
        while (argc>0) {

            // help
            if (strcmp(argv[0], "-h") == 0) {
                print_usage();
                exit(EXIT_SUCCESS);
            }

            // algorithm
            if (strcmp(argv[0], "-algo") == 0) {
                if (argc<2)
                    err("Option '-algo' expects a (string) argument.");
                algo = argv[1];
                argc-=2; argv+=2;
            }
            // ref
            else if (strcmp(argv[0], "-ref") == 0) {
                if (argc<2)
                    err("Option '-ref' expects a (filename) argument.");
                ref_filename = argv[1];
                argc-=2; argv+=2;
            }
            // flo
            else if (strcmp(argv[0], "-flo") == 0) {
                if (argc<2)
                    err("Option '-flo' expects a (filename) argument.");
                flo_filename = argv[1];
                argc-=2; argv+=2;
            }
            // ref engine
            else if (strcmp(argv[0], "-eng_ref") == 0) {
                if (argc<2)
                    err("Option '-eng_ref' expects a (string) argument.");
                eng_ref = argv[1];
                argc-=2; argv+=2;
            }
            // flo engine
            else if (strcmp(argv[0], "-eng_flo") == 0) {
                if (argc<2)
                    err("Option '-eng_flo' expects a (string) argument.");
                eng_flo = argv[1];
                argc-=2; argv+=2;
            }
            // output
            else if (strcmp(argv[0], "-output") == 0) {
                if (argc<2)
                    err("Option '-output' expects a (filename) argument.");
                output = argv[1];
                argc-=2; argv+=2;
            }
            // interpolation
            else if (strcmp(argv[0], "-interp") == 0) {
                if (argc<2)
                    err("Option '-interp' expects a (numerical) argument.");
                interp = static_cast<Resample<float>::InterpolationType>(atoi(argv[1]));
                argc-=2; argv+=2;
            }
            // add affine
            else if (strcmp(argv[0], "-add_affine") == 0) {
                if (argc<2)
                    err("Option '-add_affine' expects a (numerical) argument.");
                trans.push_back(std::make_shared<const AffineTransformation<float> >(argv[1]));
                argc-=2; argv+=2;
            }
            // add deformation
            else if (strcmp(argv[0], "-add_def") == 0) {
                if (argc<2)
                    err("Option '-add_def' expects a (numerical) argument.");
                trans.push_back(std::make_shared<const NiftiImageData3DDeformation<float> >(argv[1]));
                argc-=2; argv+=2;
            }
            // add displacement
            else if (strcmp(argv[0], "-add_disp") == 0) {
                if (argc<2)
                    err("Option '-add_disp' expects a (numerical) argument.");
                trans.push_back(std::make_shared<const NiftiImageData3DDisplacement<float> >(argv[1]));
                argc-=2; argv+=2;
            }
            // padding
            else if (strcmp(argv[0], "-pad") == 0) {
                if (argc<2)
                    err("Option '-pad' expects a (numerical) argument.");
                pad = atof(argv[1]);
                pad_set = true;
                argc-=2; argv+=2;
            }
            // direction
            else if (strcmp(argv[0], "-adj") == 0) {
                forward = false;
                argc-=1; argv+=1;
            }

            // Unknown argument
            else
                err("Unknown option '" + std::string(argv[0]) + "'.");
        }

        // Check filenames aren't blank
        if (ref_filename.size() == 0)
            err("Error: -ref required.");
        if (flo_filename.size() == 0)
            err("Error: -flo required.");

        // Get images as NiftiImages
        std::shared_ptr<const ImageData> ref = image_as_sptr(ref_filename,eng_ref);
        std::shared_ptr<const ImageData> flo = image_as_sptr(flo_filename,eng_flo);

        // Resample
        std::shared_ptr<Resample<float> > res = algo_as_sptr(algo);
        res->set_reference_image(ref);
        res->set_floating_image(flo);
        for (size_t i=0; i<trans.size(); ++i)
            res->add_transformation(trans[i]);
        res->set_interpolation_type(interp);
        if (pad_set)
            res->set_padding_value(pad);

        std::shared_ptr<ImageData> output_sptr;
        if (forward)
            output_sptr = res->forward(flo);
        else
            output_sptr = res->adjoint(ref);
        output_sptr->write(output);
    }

    // If there was an error
    catch(const std::exception &error) {
        std::cerr << "\nError encountered:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
