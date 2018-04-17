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
#include "SIRFRegActivityCorrect.h"
#include "SIRFRegImageWeightedMean.h"

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
    string wm_im2                   = examples_path + "/weighted_mean/regis_recon_gate2.nii";
    string wm_im3                   = examples_path + "/weighted_mean/regis_recon_gate3.nii";
    string wm_im4                   = examples_path + "/weighted_mean/regis_recon_gate4.nii";

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

    // ----------------------------------------------------------------------- //
    //                           Nifty aladin
    //------------------------------------------------------------------------ //
    SIRFRegNiftyAladinSym<float> NA;
    NA.set_reference_image_filename      (    reference_image_filename   );
    NA.set_floating_image_filename       (     floating_image_filename   );
    NA.set_parameter_file                (      parameter_file_aladin    );
    NA.update();
    NA.save_warped_image                 (         aladin_warped         );
    NA.save_transformation_matrix_fwrd   (             TM_fwrd           );
    NA.save_transformation_matrix_back   (             TM_back           );
    NA.save_displacement_field_fwrd_image( aladin_disp_fwrd, true,  true );
    NA.save_displacement_field_back_image( aladin_disp_back, true,  true );

    // ----------------------------------------------------------------------- //
    //                           Nifty f3d
    //------------------------------------------------------------------------ //
    SIRFRegNiftyF3dSym<float> NF;
    NF.set_reference_image_filename      (  reference_image_filename  );
    NF.set_floating_image_filename       (   floating_image_filename  );
    NF.set_parameter_file                (     parameter_file_f3d     );
    NF.set_reference_time_point          (             1              );
    NF.set_floating_time_point           (             1              );
    NF.update();
    NF.save_warped_image                 (         f3d_warped         );
    NF.save_displacement_field_fwrd_image( f3d_disp_fwrd, true,  true );
    NF.save_displacement_field_fwrd_image( f3d_disp_back, true,  true );

    // ----------------------------------------------------------------------- //
    //                           Nifty resample
    //------------------------------------------------------------------------ //
    SIRFRegNiftyResample NR;
    NR.set_reference_image_filename       ( reference_image_filename );
    NR.set_floating_image_filename        ( floating_image_filename  );
    NR.add_transformation_matrix_filename (          matrix          );
    NR.set_interpolation_type_to_cubic_spline();
    NR.update();
    NR.save_resampled_image               (       output_resample    );

    // ----------------------------------------------------------------------- //
    //                           Activity correction
    //------------------------------------------------------------------------ //
    SIRFRegActivityCorrect AC;
    AC.set_initial_activity    (        267000000         );
    AC.set_half_life           (          6586.2          );
    AC.set_input_image_filename( reference_image_filename );
    AC.set_start               (            0             );
    AC.set_stop                (           10             );
    AC.update();
    AC.save_output             (    output_activity_corr  );

    // ----------------------------------------------------------------------- //
    //                           Weighted mean
    //------------------------------------------------------------------------ //
    SIRFRegImageWeightedMean WM;
    WM.add_image         (     wm_im2, 0.2    );
    WM.add_image         (     wm_im3, 0.2    );
    WM.add_image         (     wm_im4, 0.2    );
    WM.update();
    WM.save_image_to_file(output_weighted_mean);

    // If there was an error
    } catch(const exception &error) {
        cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
