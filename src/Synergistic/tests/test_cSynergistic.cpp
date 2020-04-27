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
#include "sirf/Reg/AffineTransformation.h"
#include "sirf/Gadgetron/gadgetron_x.h"
#include "sirf/Gadgetron/gadget_lib.h"

using namespace sirf;

#define ADD_GADGET(X, T) X.push_back(std::make_shared<T>())

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

        if (argc < 3 || argc > 5) {
            std::cout << "\ntest_cSynergistic raw_mr_data nifti_filename [mr_recon_h5_filename]\n";
            return EXIT_SUCCESS;
        }

        const std::string raw_mr_filename = argv[1];
        const std::string nifti_filename = argv[2];
        std::string mr_recon_h5_filename = "";
        if (argc > 3)
            mr_recon_h5_filename = argv[3];

        // Test STIR -> Nifti
        {
            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Starting STIRImageData->NiftiImageData test...\n";
            std::cout << "//------------------------------------------------------------------------ //\n";

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

            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Finished STIRImageData->NiftiImageData test.\n";
            std::cout << "//------------------------------------------------------------------------ //\n";
        }

        // Test Gadgetron -> Nifti
        if (!mr_recon_h5_filename.empty()) {

            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Starting GadgetronImageData->NiftiImageData test...\n";
            std::cout << "//------------------------------------------------------------------------ //\n";

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

            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Finished GadgetronImageData->NiftiImageData test.\n";
            std::cout << "//------------------------------------------------------------------------ //\n";
        }
        else {
            std::cout << "\nNot performing GadgetronImagesVector to NiftiImageData conversion (as .h5 file was not provided).\n";
        }   

        // Test complex resample (need to do a reconstruction first)
        {
            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Starting complex resampler test...\n";
            std::cout << "//------------------------------------------------------------------------ //\n";

            AcquisitionsFile raw_mr(raw_mr_filename);

            std::vector<gadgetron::shared_ptr<Gadget> > gadgets;
            ADD_GADGET(gadgets, NoiseAdjustGadget);
            ADD_GADGET(gadgets, AsymmetricEchoAdjustROGadget);
            ADD_GADGET(gadgets, RemoveROOversamplingGadget);
            ADD_GADGET(gadgets, AcquisitionAccumulateTriggerGadget);
            ADD_GADGET(gadgets, BucketToBufferGadget);
            ADD_GADGET(gadgets, SimpleReconGadget);
            ADD_GADGET(gadgets, ImageArraySplitGadget);

            ImagesReconstructor recon;
            for (unsigned i=0; i<gadgets.size(); ++i)
                recon.add_gadget("gadget_" + std::to_string(i), gadgets[i]);
            recon.process(raw_mr);
            std::shared_ptr<GadgetronImageData> ismrmrd_im_sptr = recon.get_output();

            if (!ismrmrd_im_sptr->is_complex())
                throw std::runtime_error("Expected output of reconstruction to be complex");

            // Complex component may be empty, so use imag = real / 2
            if ((*ismrmrd_im_sptr->begin()).get_typeID() != NumberType::CXFLOAT)
                throw std::runtime_error("Expected output of reconstruction to be complex float");
            for (auto &iter = ismrmrd_im_sptr->begin(); iter != ismrmrd_im_sptr->end(); ++iter) {
                complex_float_t cmplx_flt(*iter);
                cmplx_flt.imag(cmplx_flt.real() / 2.f);
                *iter = NumRef((void *)&cmplx_flt, NumberType::CXFLOAT);
            }

            // Convert the complex image to two niftis
            std::shared_ptr<NiftiImageData<float> > real_sptr, imag_sptr;
            NiftiImageData<float>::construct_NiftiImageData_from_complex_im(real_sptr,imag_sptr,ismrmrd_im_sptr);

            real_sptr->write("results/real");
            imag_sptr->write("results/imag");

            // Create affine transformation
            std::shared_ptr<AffineTransformation<float> > tm_sptr =
                    std::make_shared<AffineTransformation<float> >();
            (*tm_sptr)[0][3] = 2.f;

            // Resample the complex data
            NiftyResample<float> res_complex;
            res_complex.set_reference_image(ismrmrd_im_sptr);
            res_complex.set_floating_image(ismrmrd_im_sptr);
            res_complex.set_interpolation_type_to_linear();
            res_complex.add_transformation(tm_sptr);
            std::shared_ptr<ImageData> forward_cplx_sptr = res_complex.forward(ismrmrd_im_sptr);
            std::shared_ptr<ImageData> adjoint_cplx_sptr = res_complex.adjoint(ismrmrd_im_sptr);

            // Get the output
            std::shared_ptr<NiftiImageData<float> > forward_cplx_real_sptr, forward_cplx_imag_sptr, adjoint_cplx_real_sptr, adjoint_cplx_imag_sptr;
            NiftiImageData<float>::construct_NiftiImageData_from_complex_im(forward_cplx_real_sptr,forward_cplx_imag_sptr,forward_cplx_sptr);
            NiftiImageData<float>::construct_NiftiImageData_from_complex_im(adjoint_cplx_real_sptr,adjoint_cplx_imag_sptr,adjoint_cplx_sptr);

            forward_cplx_real_sptr->write("results/forward_cplx_real");
            forward_cplx_imag_sptr->write("results/forward_cplx_imag");
            adjoint_cplx_real_sptr->write("results/adjoint_cplx_real");
            adjoint_cplx_imag_sptr->write("results/adjoint_cplx_imag");

            // Now resample each of the components individually
            NiftyResample<float> res_real;
            res_real.set_reference_image(real_sptr);
            res_real.set_floating_image(real_sptr);
            res_real.set_interpolation_type_to_linear();
            res_real.add_transformation(tm_sptr);
            std::shared_ptr<NiftiImageData<float> > forward_real_sptr =
                    std::dynamic_pointer_cast<NiftiImageData<float> >(res_real.forward(real_sptr));
            std::shared_ptr<NiftiImageData<float> > adjoint_real_sptr =
                    std::dynamic_pointer_cast<NiftiImageData<float> >(res_real.adjoint(real_sptr));

            NiftyResample<float> res_imag;
            res_imag.set_reference_image(imag_sptr);
            res_imag.set_floating_image(imag_sptr);
            res_imag.set_interpolation_type_to_linear();
            res_imag.add_transformation(tm_sptr);
            std::shared_ptr<NiftiImageData<float> > forward_imag_sptr =
                    std::dynamic_pointer_cast<NiftiImageData<float> >(res_imag.forward(imag_sptr));
            std::shared_ptr<NiftiImageData<float> > adjoint_imag_sptr =
                    std::dynamic_pointer_cast<NiftiImageData<float> >(res_imag.adjoint(imag_sptr));

            // Compare that the real and imaginary parts match regardless
            // of whether they were resampled separately or together.
            if (*forward_real_sptr != *forward_cplx_real_sptr
                    || *forward_imag_sptr != *forward_cplx_imag_sptr)
                throw std::runtime_error("NiftyResample forward failed for complex data");
            if (*adjoint_real_sptr != *adjoint_cplx_real_sptr
                    || *adjoint_imag_sptr != *adjoint_cplx_imag_sptr)
                throw std::runtime_error("NiftyResample adjoint failed for complex data");

            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Finished complex resampler test.\n";
            std::cout << "//------------------------------------------------------------------------ //\n";
        }

        // Test MR reorient
        {
            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Starting GadgetronImageData reorient test...\n";
            std::cout << "//------------------------------------------------------------------------ //\n";

            // Read ISMRMRD image
            std::shared_ptr<GadgetronImagesVector> G1_sptr = std::make_shared<GadgetronImagesVector>();
            G1_sptr->read(mr_recon_h5_filename);

            // Convert ISMRMRD image to nifti
            std::shared_ptr<NiftiImageData<float> > G1_nii_sptr =
                    std::make_shared<NiftiImageData<float> >(*G1_sptr);
            std::shared_ptr<nifti_image> raw_nii_sptr = G1_nii_sptr->get_raw_nifti_sptr();

            // Affine transformation as translation by integer num voxels (so no interpolation)
            std::array<float,3> trans = {  G1_sptr->get_geom_info_sptr()->get_spacing()[0] * 2.f,
                                           G1_sptr->get_geom_info_sptr()->get_spacing()[1] * 5.f,
                                           G1_sptr->get_geom_info_sptr()->get_spacing()[2] * 3.f };

            // Create matching transformation matrix
            std::shared_ptr<AffineTransformation<float> > trans_sptr =
                    std::make_shared<AffineTransformation<float> >();
            for (unsigned i=0; i<3; ++i)
                (*trans_sptr)[i][3] = trans[i];

            // Shift nifti
            for (unsigned i=0; i<3; ++i)
                raw_nii_sptr->qto_xyz.m[i][3] += trans[i];
            raw_nii_sptr->qto_ijk =
                    nifti_mat44_inverse(raw_nii_sptr->qto_xyz);
            nifti_mat44_to_quatern( raw_nii_sptr->qto_xyz,
                                    &raw_nii_sptr->quatern_b,
                                    &raw_nii_sptr->quatern_c,
                                    &raw_nii_sptr->quatern_d,
                                    &raw_nii_sptr->qoffset_x,
                                    &raw_nii_sptr->qoffset_y,
                                    &raw_nii_sptr->qoffset_z,
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    &raw_nii_sptr->qfac );
            raw_nii_sptr->pixdim[0]=raw_nii_sptr->qfac;
            // Check everything is ok
            reg_checkAndCorrectDimension(raw_nii_sptr.get());
            // Re-set up geom info
            G1_nii_sptr->set_up_geom_info();

            // Reorient gadgetron image with modified G1_nii_sptr's geom info
            std::shared_ptr<GadgetronImagesVector> G2_sptr = G1_sptr->clone();
            G2_sptr->reorient(*G1_nii_sptr->get_geom_info_sptr());

            std::cout << "\n original:\n";
            G1_sptr->get_geom_info_sptr()->print_info();
            std::cout << "\n resampled:\n";
            G2_sptr->get_geom_info_sptr()->print_info();

            // Now resampled G2 back to G1 using inverse TM, should be the same
            NiftyResample<float> res;
            res.set_reference_image(G1_sptr);
            res.set_floating_image(G2_sptr);
            res.set_padding_value(0.f);
            res.set_interpolation_type_to_linear();
            res.add_transformation(std::make_shared<const AffineTransformation<float> >(trans_sptr->get_inverse()));
            std::shared_ptr<ImageData> resampled_G2_sptr = res.forward(G2_sptr);

            std::cout << "\n reoriented back to original space:\n";
            resampled_G2_sptr->get_geom_info_sptr()->print_info();

            if (NiftiImageData<float>(*G1_sptr) != NiftiImageData<float>(*resampled_G2_sptr))
                throw std::runtime_error("GadgetronImagesVector::reorient test failed");


            std::cout << "// ----------------------------------------------------------------------- //\n";
            std::cout << "//                  Finished GadgetronImageData reorient test.\n";
            std::cout << "//------------------------------------------------------------------------ //\n";
        }

    // Error handling
    } catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
