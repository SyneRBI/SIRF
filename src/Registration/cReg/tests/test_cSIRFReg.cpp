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
\brief Test for creation of a weighted-mean image

\author Richard Brown
\author CCP PETMR
*/

#include <iostream>
#include "SIRFRegNiftyAladinSym.h"
#include "SIRFRegNiftyF3dSym.h"
#include "SIRFRegNiftyResample.h"
#include "SIRFImageData.h"
#include "SIRFRegImageWeightedMean.h"
#include "stir_data_containers.h"

using namespace std;

int main(int argc, char* argv[])
{

    try {

    // Get current working directory
    boost::filesystem::path path( boost::filesystem::initial_path<boost::filesystem::path>() );
    path = boost::filesystem::system_complete( boost::filesystem::path( argv[0] ) );
    path = path.remove_filename();

    // Paths
    string SIRF_PATH     = getenv("SIRF_PATH");
    string examples_path = SIRF_PATH + "/data/examples/Registration";
    string output_path   = path.string() + "/results/";

    // Input filenames
    string reference_image_filename = examples_path + "/mouseFixed.nii.gz";
    string floating_image_filename  = examples_path + "/mouseMoving.nii.gz";
    string parameter_file_aladin    = examples_path + "/paramFiles/aladin.par";
    string parameter_file_f3d       = examples_path + "/paramFiles/f3d.par";
    string matrix                   = examples_path + "/transformation_matrix.txt";
    string stir_nifti               = examples_path + "/nifti_created_by_stir.nii";

    // Output filenames
    string aladin_warped            = output_path   + "cplusplus_aladin_warped";
    string f3d_warped               = output_path   + "cplusplus_f3d_warped";
    string TM_fwrd                  = output_path   + "cplusplus_TM_fwrd.txt";
    string TM_back                  = output_path   + "cplusplus_TM_back.txt";
    string aladin_disp_fwrd         = output_path   + "cplusplus_aladin_disp_fwrd";
    string aladin_disp_back         = output_path   + "cplusplus_aladin_disp_back";
    string f3d_disp_fwrd            = output_path   + "cplusplus_f3d_disp_fwrd";
    string f3d_disp_back            = output_path   + "cplusplus_f3d_disp_back";

    string output_resample          = output_path   + "cplusplus_resample";
    string output_activity_corr     = output_path   + "cplusplus_activity_corr";
    string output_weighted_mean     = output_path   + "cplusplus_weighted_mean";

    string output_stir_nifti        = output_path   + "cplusplus_stir_nifti.nii";

    SIRFImageData reference( reference_image_filename );
    SIRFImageData floating (  floating_image_filename );
    SIRFImageData nifti    (        stir_nifti        );



    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty aladin test...                          //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegNiftyAladinSym<float> NA;
    NA.set_reference_image               (            reference          );
    NA.set_floating_image                (            floating           );
    NA.set_parameter_file                (      parameter_file_aladin    );
    NA.update();
    NA.save_warped_image                 (         aladin_warped         );
    NA.save_transformation_matrix_fwrd   (             TM_fwrd           );
    NA.save_transformation_matrix_back   (             TM_back           );
    NA.save_displacement_field_fwrd_image( aladin_disp_fwrd, true        );
    NA.save_displacement_field_back_image( aladin_disp_back, true        );
    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished Nifty aladin test.                            //\n";
    cout << "//------------------------------------------------------------------------ //\n";



    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty f3d test...                             //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegNiftyF3dSym<float> NF;
    NF.set_reference_image               (         reference          );
    NF.set_floating_image                (          floating          );
    NF.set_parameter_file                (     parameter_file_f3d     );
    NF.set_reference_time_point          (             1              );
    NF.set_floating_time_point           (             1              );
    NF.update();
    NF.save_warped_image                 (         f3d_warped         );
    NF.save_displacement_field_fwrd_image( f3d_disp_fwrd, true        );
    NF.save_displacement_field_fwrd_image( f3d_disp_back, true        );
    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished Nifty f3d test.                               //\n";
    cout << "//------------------------------------------------------------------------ //\n";




    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty resample test...                        //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegNiftyResample NR;
    NR.set_reference_image                (         reference        );
    NR.set_floating_image                 (         floating         );
    NR.set_transformation_matrix          (          matrix          );
    NR.set_interpolation_type_to_cubic_spline();
    NR.update();
    NR.save_resampled_image               (       output_resample    );
    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished Nifty resample test.                          //\n";
    cout << "//------------------------------------------------------------------------ //\n";




    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting weighted mean test...                         //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegImageWeightedMean WM;
    WM.add_image         (     nifti, 0.2F    );
    WM.add_image         (     nifti, 0.2F    );
    WM.add_image         (     nifti, 0.2F    );
    WM.update();
    WM.save_image_to_file(output_weighted_mean);
    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished weighted mean test.                           //\n";
    cout << "//------------------------------------------------------------------------ //\n";





    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting PET SIRFImageData test...                     //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    // Open stir image and convert to SIRFImageData
    sirf::PETImageData pet_image_data(stir_nifti);
    SIRFImageData image_data_from_stir(pet_image_data);
    // Compare to nifti IO (if they don't match, you'll see a message but don't throw an error for now)
    SIRFImageData image_data_from_nifti(stir_nifti);
    SIRFRegMisc::do_nift_image_match(image_data_from_stir, image_data_from_nifti);
    // Print info
    std::vector<SIRFImageData> ims;
    ims.push_back(image_data_from_stir);
    ims.push_back(image_data_from_nifti);
    SIRFRegMisc::dump_nifti_info(ims);
    // Save the one opened by stir
    image_data_from_stir.save_to_file(output_stir_nifti);
    // Now clone the converted and fill with 1's
    sirf::PETImageData cloned = pet_image_data;
    cloned.fill(1.F);
    // Fill the cloned image with data from converted
    image_data_from_stir.copy_data_to(cloned);
    // Compare
    if (fabs(cloned.data()[0][0][0] - pet_image_data.data()[0][0][0]) > 1.e-7F) {
        throw std::runtime_error("Image was not filled from nifti.");
    }
    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished PET SIRFImageData test.                       //\n";
    cout << "//------------------------------------------------------------------------ //\n";



    // Error handling
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
