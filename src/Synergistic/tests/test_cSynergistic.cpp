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
\ingroup Synergistic
\brief Synergistic tests

\author Richard Brown
\author CCP PETMR
*/

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Reg/NiftyResample.h"

using namespace sirf;

static void create_stir_output_file_format(const std::string &path)
{
    std::ofstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Unable to write stir output file format.");
    file << "OutputFileFormat Parameters:=\n";
    file << "output file format type := ITK\n";
    file << "ITK Output File Format Parameters:=\n";
    file << "number format := float\n";
    file << "number_of_bytes_per_pixel:=4\n";
    file << "default extension:=.nii\n";
    file << "End ITK Output File Format Parameters:=\n";
    file << "End:=\n";
    file.close();
}

int main(int argc, char* argv[])
{
    try {

        if (argc < 2) {
            std::cout << "\ntest_cSynergistic nifti_filename [mr_recon_h5_filename]\n";
            return EXIT_SUCCESS;
        }

        const std::string nifti_filename = argv[1];
        std::string mr_recon_h5_filename = "";
        if (argc > 2)
            mr_recon_h5_filename = argv[2];

        // Test STIR -> Nifti
        {
            std::cout << "\nPerforming STIRImageData to NiftiImageData conversion...\n";

            // Load the image as a NiftiImageData3D
            NiftiImageData3D<float> image_nifti(nifti_filename);

            // Read as STIRImageData, convert NiftiImageData3D and save to file
            STIRImageData image_stir(nifti_filename);
            NiftiImageData3D<float> image_nifti_from_stir(image_stir);
            image_nifti_from_stir.write("results/stir_to_nifti.nii",image_nifti.get_original_datatype());

            // Compare the two
            if (image_nifti != image_nifti_from_stir)
                throw std::runtime_error("Conversion from STIR to Nifti failed");

            // Also save the STIRImageData to file (might be useful visual for comparison)
            create_stir_output_file_format("results/stir_output_file_format_nifti.par");
            image_stir.write("results/stir.nii","results/stir_output_file_format_nifti.par");
        }

        // Test Gadgetron -> Nifti
        if (!mr_recon_h5_filename.empty()) {

            std::cout << "\nPerforming GadgetronImagesVector to NiftiImageData conversion...\n";

            // Read ISMRMRD image
            GadgetronImagesVector ismrmrd_im;
            ismrmrd_im.read(mr_recon_h5_filename);

            // Convert ISMRMRD image to nifti
            NiftiImageData<float> nifti_from_ismrmrd(ismrmrd_im);

            // Read vendor-reconstructed image
            NiftiImageData<float> dicom_im(nifti_filename);

            // Standardise to remove scaling problems
            dicom_im.standardise();
            nifti_from_ismrmrd.standardise();

            // Compare the two. Since the images are being reconstructed independently, there is no
            // guarantee they will perfectly match. So we need an data-dependent acceptance threshold.
            if (!NiftiImageData<float>::are_equal_to_given_accuracy(dicom_im, nifti_from_ismrmrd, 165.f))
                throw std::runtime_error("Conversion from ISMRMRD to Nifti failed");
        }
        else {
            std::cout << "\nNot performing GadgetronImagesVector to NiftiImageData conversion (as .h5 file was not provided).\n";
        }

    // Error handling
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
