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
#include <memory>

using namespace std;

int main(int argc, char* argv[])
{

    try {

    // Paths
    string SIRF_PATH;
    if (argc==1)
        SIRF_PATH = getenv("SIRF_PATH");
    else
        SIRF_PATH = argv[1];
    string examples_path = SIRF_PATH + "/data/examples/Registration";
    string output_path   = "results/";

    // Input filenames
    string ref_aladin_filename      = examples_path + "/test.nii.gz";
    string flo_aladin_filename      = examples_path + "/test2.nii.gz";
    string ref_f3d_filename         = examples_path + "/mouseFixed.nii.gz";
    string flo_f3d_filename         = examples_path + "/mouseMoving.nii.gz";
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
    string rigid_resample           = output_path   + "cplusplus_rigid_resample";
    string nonrigid_resample_def    = output_path   + "cplusplus_nonrigid_resample_def";
    string nonrigid_resample_disp   = output_path   + "cplusplus_nonrigid_resample_disp";
    string output_weighted_mean     = output_path   + "cplusplus_weighted_mean";
    string output_weighted_mean_def = output_path   + "cplusplus_weighted_mean_def";

    string output_stir_nifti        = output_path   + "cplusplus_stir_nifti.nii";

    SIRFImageData ref_aladin( ref_aladin_filename );
    SIRFImageData flo_aladin( flo_aladin_filename );
    SIRFImageData ref_f3d   (   ref_f3d_filename  );
    SIRFImageData flo_f3d   (   flo_f3d_filename  );
    SIRFImageData nifti     (      stir_nifti     );

    float required_percentage_accuracy = 1.F;


    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty aladin test...                          //\n";
    cout << "//------------------------------------------------------------------------ //\n";

    SIRFRegNiftyAladinSym<float> NA;
    NA.set_reference_image               (           ref_aladin          );
    NA.set_floating_image                (           flo_aladin          );
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
    NF.set_reference_image               (          ref_f3d           );
    NF.set_floating_image                (          flo_f3d           );
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
    cout << "//                  Starting transformations test...                       //\n";
    cout << "//------------------------------------------------------------------------ //\n";

    cout << "\nchecking if affine to deformation is ok...\n";
    SIRFRegTransformationAffine trans_affine(NA.get_transformation_matrix_fwrd());
    if (!SIRFRegMisc::do_nifti_images_match(trans_affine.get_as_deformation_field(ref_aladin), NA.get_deformation_field_fwrd(), required_percentage_accuracy))
        throw runtime_error("Affine as deformation does not match deformation.");

    cout << "\nchecking if displacement to deformation is ok...\n";
    SIRFRegTransformationDisplacement trans_disp(NA.get_displacement_field_fwrd());
    if (!SIRFRegMisc::do_nifti_images_match(trans_disp.get_as_deformation_field(ref_aladin), NA.get_deformation_field_fwrd(), required_percentage_accuracy))
        throw runtime_error("Disp as deformation does not match deformation.");

    cout << "\nchecking if deformation to deformation is ok...\n";
    SIRFRegTransformationDeformation  trans_def(NA.get_deformation_field_fwrd());
    if (!SIRFRegMisc::do_nifti_images_match(trans_def.get_as_deformation_field(ref_aladin), NA.get_deformation_field_fwrd(), required_percentage_accuracy))
        throw runtime_error("Def as deformation does not match deformation.");

    cout << "\njoining two identity matrices and the disp field. convert to deformation, should be the same as the deformation created by the same registration...\n";
    // Check that the composition of two identity matrices and the displacement field equals the deformation field image.
    mat44 identity = SIRFRegMisc::get_identity_matrix();
    std::vector<std::shared_ptr<SIRFRegTransformation> > vec;
    vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegTransformationAffine(identity)));
    vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegTransformationAffine(identity)));
    vec.push_back(std::shared_ptr<SIRFRegTransformation>(new SIRFRegTransformationDisplacement(trans_disp)));
    SIRFRegTransformationDeformation trans_composed;
    SIRFRegMisc::compose_transformations_into_single_deformation(trans_composed, vec, ref_aladin);
    if (!SIRFRegMisc::do_nifti_images_match(trans_composed.get_as_deformation_field(ref_aladin), NA.get_deformation_field_fwrd(), required_percentage_accuracy))
        throw runtime_error("identity(identity(displacement)) as deformation does not match deformation.");

    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished transformations test.                         //\n";
    cout << "//------------------------------------------------------------------------ //\n";




    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty resample rigid test...                  //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegNiftyResample NRA;
    NRA.set_reference_image                (        ref_aladin        );
    NRA.set_floating_image                 (        flo_aladin        );
    NRA.add_transformation_affine          (          identity        );
    NRA.add_transformation_affine          (          TM_fwrd         );
    NRA.set_interpolation_type_to_cubic_spline();
    NRA.update();
    NRA.save_resampled_image               (        rigid_resample    );

    if (!SIRFRegMisc::do_nifti_images_match(NA.get_output(),NRA.get_output(), required_percentage_accuracy))
        throw runtime_error("Resampled image (rigid) does not match the one resampled with aladin.");

    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished Nifty resample rigid test.                    //\n";
    cout << "//------------------------------------------------------------------------ //\n";




    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty resample non-rigid deformation test...  //\n";
    cout << "//------------------------------------------------------------------------ //\n";

    SIRFRegNiftyResample NRF2;
    NRF2.set_reference_image                (               ref_f3d            );
    NRF2.set_floating_image                 (               flo_f3d            );
    NRF2.set_interpolation_type_to_cubic_spline();
    NRF2.add_transformation_def             (  NF.get_deformation_field_fwrd() );
    NRF2.update();
    NRF2.save_resampled_image               (       nonrigid_resample_def      );

    if (!SIRFRegMisc::do_nifti_images_match(NF.get_output(),NRF2.get_output(), required_percentage_accuracy))
        throw runtime_error("Resampled image (non-rigid) does not match the one resampled with f3d.");

    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished Nifty resample non-rigid deformation test.    //\n";
    cout << "//------------------------------------------------------------------------ //\n";




    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting Nifty resample non-rigid displacement test... //\n";
    cout << "//------------------------------------------------------------------------ //\n";

    SIRFRegNiftyResample NRF1;
    NRF1.set_reference_image                (               ref_f3d            );
    NRF1.set_floating_image                 (               flo_f3d            );
    NRF1.set_interpolation_type_to_cubic_spline();
    NRF1.add_transformation_disp            ( NF.get_displacement_field_fwrd() );
    NRF1.update();
    NRF1.save_resampled_image               (      nonrigid_resample_disp      );

    if (!SIRFRegMisc::do_nifti_images_match(NF.get_output(),NRF1.get_output(), required_percentage_accuracy))
        throw runtime_error("Resampled image (non-rigid) does not match the one resampled with f3d.");

    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished Nifty resample non-rigid displacement test.   //\n";
    cout << "//------------------------------------------------------------------------ //\n";





    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting weighted mean test...                         //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegImageWeightedMean<SIRFImageData> WM;
    SIRFImageData im1(stir_nifti);
    SIRFImageData im2(stir_nifti);
    SIRFImageData im3(stir_nifti);
    SIRFImageData im4(stir_nifti);
    im1.fill(1.F);
    im2.fill(4.F);
    im3.fill(7.F);
    im4.fill(6.F);

    WM.add_image( im1, 2.F );
    WM.add_image( im2, 4.F );
    WM.add_image( im3, 3.F );
    WM.add_image( im4, 1.F );
    WM.update();
    WM.save_image_to_file(output_weighted_mean);

    // Answer should be 4.5, so compare it to that!
    SIRFImageData res(stir_nifti);
    res.fill(4.5F);

    SIRFImageData dt = WM.get_output();
    //SIRFImageDataDeformation def3 = dynamic_cast<SIRFImageDataDeformation&>(dt);

    if (!SIRFRegMisc::do_nifti_images_match(WM.get_output(), res, required_percentage_accuracy))
        throw runtime_error("Weighted mean does not match the original.");

    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished weighted mean test.                           //\n";
    cout << "//------------------------------------------------------------------------ //\n";





    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting weighted mean deformation test...             //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    SIRFRegImageWeightedMean<SIRFImageDataDeformation> WM_def;
    WM_def.add_image         (     NF.get_deformation_field_fwrd(), 0.2F    );
    WM_def.add_image         (     NF.get_deformation_field_fwrd(), 0.2F    );
    WM_def.add_image         (     NF.get_deformation_field_fwrd(), 0.2F    );
    WM_def.update();
    WM_def.save_image_to_file(output_weighted_mean_def);
    SIRFImageDataDeformation def = WM_def.get_output();

    if (!SIRFRegMisc::do_nifti_images_match(NF.get_deformation_field_fwrd() ,WM_def.get_output(), required_percentage_accuracy))
        throw runtime_error("Weighted mean does not match the original.");

    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Finished weighted mean deformation test.               //\n";
    cout << "//------------------------------------------------------------------------ //\n";





    cout << "// ----------------------------------------------------------------------- //\n";
    cout << "//                  Starting PET SIRFImageData test...                     //\n";
    cout << "//------------------------------------------------------------------------ //\n";
    // Open stir image and convert to SIRFImageData
    sirf::PETImageData pet_image_data(stir_nifti);
    SIRFImageData image_data_from_stir(pet_image_data);
    // Compare to nifti IO (if they don't match, you'll see a message but don't throw an error for now)
    SIRFImageData image_data_from_nifti(stir_nifti);
    SIRFRegMisc::do_nifti_images_match(image_data_from_stir, image_data_from_nifti, required_percentage_accuracy);
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
