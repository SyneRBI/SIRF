/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2019 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2019 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
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
\brief MR related C++ tests

\author Johannes Mayer
\author SyneRBI
*/

#include <iostream>
#include <cstdlib>
#include <numeric>
#include <vector>


#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"

#include "sirf/Gadgetron/encoding.h"
#include "sirf/Gadgetron/chain_lib.h"

#include "mrtest_auxiliary_funs.h"

using namespace sirf;


bool test_TrajectoryPreparation_constructors( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::CartesianTrajectoryPrep cart_tp;
        sirf::GRPETrajectoryPrep rpe_tp;

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}


bool test_GRPETrajectoryPrep_set_trajectory(const AcquisitionsVector av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        sirf::GRPETrajectoryPrep rpe_tp;

        AcquisitionsVector av_temp(av);

        rpe_tp.set_trajectory(av_temp);
        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}


bool test_get_kspace_order(const MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;


        auto kspace_sorting_slice = av.get_kspace_order();

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_get_subset(const MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        std::vector<int> subset_idx;
        for(int i=0; i<av.number()/10; ++i)
            subset_idx.push_back(i);

        sirf::AcquisitionsVector subset;
        av.get_subset(subset, subset_idx);

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_CoilSensitivitiesVector_calculate( MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        CoilSensitivitiesVector csv;
        csv.set_csm_smoothness(50);
        csv.calculate(av);

        std::cout << "We have " << csv.items() << " coilmaps" << std::endl;

        for(int i=0; i<csv.items(); ++i)
        {
            gadgetron::shared_ptr<ImageWrap> sptr_iw = csv.sptr_image_wrap(i);

            std::stringstream fname_out;
            fname_out << "output_" << __FUNCTION__ << "_" << i;

            sirf::write_cfimage_to_raw(fname_out.str(), *sptr_iw);
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_CoilSensitivitiesVector_get_csm_as_cfimage(MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        CoilSensitivitiesVector csv;
        csv.calculate(av);

        std::cout << "We have " << csv.items() << " coilmaps" << std::endl;

        for(int i=0; i<csv.items(); ++i)
        {
            CFImage img = csv.get_csm_as_cfimage(i);

            std::stringstream fname_out;
            fname_out << "output_" << __FUNCTION__ << "_" << i;

            sirf::write_cfimage_to_raw(fname_out.str(), img);
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_number_of_encodings(gadgetron::shared_ptr<MRAcquisitionData> sptr_ad)
{
    MRAcquisitionData& ad = *sptr_ad;
    gadgetron::shared_ptr<CoilSensitivitiesVector> sptr_csv
        (new CoilSensitivitiesVector);
    CoilSensitivitiesVector& csv = *sptr_csv;
    csv.calculate(ad);

    gadgetron::shared_ptr<GadgetronImageData> sptr_id;
    gadgetron::shared_ptr<GadgetronImageData> sptr;

    std::cout << "################ Directorly BEFORE calling undersampled(): # encodings: " << ad.acquisitions_info().get_IsmrmrdHeader().encoding.size() << " encodings now ####################### " << std::endl;
    ad.undersampled();

    if (ad.undersampled()) {
        std::cout << "found undersampled data\n";
        SimpleGRAPPAReconstructionProcessor recon;
        std::cout << "reconstructing...\n";
        recon.process(ad);
        sptr_id = recon.get_output();
        sptr_id = sptr_id->clone("GADGETRON_DataRole", "image");
    }
    else {
        std::cout << "found fully sampled data\n";
        SimpleReconstructionProcessor recon;
        std::cout << "reconstructing...\n";
        recon.process(ad);
        sptr_id = recon.get_output();
    }
}

bool test_acq_mod_norm(gadgetron::shared_ptr<MRAcquisitionData> sptr_ad)
{
    MRAcquisitionData& ad = *sptr_ad;
    int acq_dim[10];
    ad.get_acquisitions_dimensions((size_t)acq_dim);
    std::cout << "acquisitions dimensions: "
              << acq_dim[0] << " by "
              << acq_dim[1] << " by "
              << acq_dim[2] << '\n';

    gadgetron::shared_ptr<CoilSensitivitiesVector> sptr_csv
        (new CoilSensitivitiesVector);
    CoilSensitivitiesVector& csv = *sptr_csv;
    csv.calculate(ad);

    gadgetron::shared_ptr<GadgetronImageData> sptr_id;
    gadgetron::shared_ptr<GadgetronImageData> sptr;
    if (ad.undersampled()) {
        std::cout << "found undersampled data\n";
        SimpleGRAPPAReconstructionProcessor recon;
        std::cout << "reconstructing...\n";
        recon.process(ad);
        sptr_id = recon.get_output();
        sptr_id = sptr_id->clone("GADGETRON_DataRole", "image");
    }
    else {
        std::cout << "found fully sampled data\n";
        SimpleReconstructionProcessor recon;
        std::cout << "reconstructing...\n";
        recon.process(ad);
        sptr_id = recon.get_output();
    }

    GadgetronImageData& image = *sptr_id;
    int img_dim[10];
    int ni = image.number();
    image.get_image_dimensions(0, img_dim);
    std::cout << ni << " images reconstructed\n";
    std::cout << "image dimensions: "
              << img_dim[0] << " by "
              << img_dim[1] << " by "
              << img_dim[2] << " by "
              << img_dim[3] << '\n';

    MRAcquisitionModel am;
    am.set_up(sptr_ad, sptr_id);
    am.setCSMs(sptr_csv);
    std::cout << "forward projectiong...\n";
    gadgetron::shared_ptr<MRAcquisitionData> sptr_sd = am.fwd(image);
    MRAcquisitionData& sd = *sptr_sd;
    sd.get_acquisitions_dimensions((size_t)acq_dim);
    std::cout << "simulated acquisitions dimensions: "
              << acq_dim[0] << " by "
              << acq_dim[1] << " by "
              << acq_dim[2] << '\n';

    float im_norm = image.norm();
    float sd_norm = sd.norm();
    float am_norm = am.norm();
    std::cout << "\nchecking the acquisition model norm:\n";
    std::cout << "acquisition model norm: |A| = " << am_norm << '\n';
    std::cout << "image data x norm: |x| = " << im_norm << '\n';
    std::cout << "simulated acquisition data norm: |A(x)| = " << sd_norm << '\n';
    std::cout << "checking that |A(x)| <= |A||x|: ";
    float bound = am_norm*im_norm;
    bool ok = (sd_norm <= bound);
    if (ok)
        std::cout << sd_norm << " <= " << bound << " ok!\n";
    else
        std::cout << sd_norm << " > " << bound << " failure!\n";

    return ok;
}

bool test_bwd(MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::GadgetronImagesVector img_vec;
        sirf::MRAcquisitionModel acquis_model;

        sirf::CoilSensitivitiesVector csm;
        csm.calculate(av);

        auto sptr_encoder = std::make_shared<sirf::CartesianFourierEncoding>(sirf::CartesianFourierEncoding());
        acquis_model.set_encoder(sptr_encoder);

        acquis_model.bwd(img_vec, csm, av);

        for(int i=0; i<img_vec.items(); ++i)
        {
            std::stringstream fname_output;
            fname_output << "output_" << __FUNCTION__ << "_image_" << i;
            write_cfimage_to_raw(fname_output.str(), img_vec.image_wrap(i));
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_get_rpe_trajectory(AcquisitionsVector av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        sirf::GRPETrajectoryPrep rpe_tp;

        rpe_tp.set_trajectory(av);

        std::stringstream fname_output;
        fname_output << "output_" << __FUNCTION__ << ".h5";
        av.write(fname_output.str());

        std::cout << "wrote file apparently" << std::endl;

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}


bool test_rpe_csm(MRAcquisitionData& av)
{
    try
    {
       std::cout << "Running test " << __FUNCTION__ << std::endl;


       sirf::GRPETrajectoryPrep rpe_tp;
       rpe_tp.set_trajectory(av);

       sirf::CoilSensitivitiesVector csm;
       csm.set_csm_smoothness(50);
       csm.calculate(av);

       for(int i=0; i<csm.items(); ++i)
       {
           std::stringstream fname_output;
           fname_output << "output_" << __FUNCTION__ << "_image_" << i;
           write_cfimage_to_raw(fname_output.str(), csm.image_wrap(i));
       }

       return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }




}


bool test_rpe_bwd(MRAcquisitionData& av)
{
    try
    {
       std::cout << "Running test " << __FUNCTION__ << std::endl;


       sirf::GRPETrajectoryPrep rpe_tp;
       rpe_tp.set_trajectory(av);


       sirf::GadgetronImagesVector img_vec;
       img_vec.set_meta_data(av.acquisitions_info());

       if(!av.sorted())
           av.sort();

       auto sort_idx = av.get_kspace_order();

       auto sptr_enc = std::make_shared< RPEFourierEncoding > (RPEFourierEncoding());

       for(int i=0; i<sort_idx.size(); ++i)
       {
           sirf::AcquisitionsVector subset;
           av.get_subset(subset, sort_idx[i]);

           CFImage img;
           sptr_enc->backward(&img, subset);

           void* vptr_img = new CFImage(img);// god help me I don't trust this!
           ImageWrap iw(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, vptr_img);

           img_vec.append(iw);

       }

       for(int i=0; i<img_vec.items(); ++i)
       {
           std::stringstream fname_output;
           fname_output << "output_" << __FUNCTION__ << "_image_" << i;
           write_cfimage_to_raw(fname_output.str(), img_vec.image_wrap(i));
       }

       return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}


bool test_rpe_fwd(MRAcquisitionData& av)
{
    try
    {
       std::cout << "Running test " << __FUNCTION__ << std::endl;

       sirf::GRPETrajectoryPrep rpe_tp;
       rpe_tp.set_trajectory(av);

       sirf::GadgetronImagesVector img_vec;
       img_vec.set_meta_data(av.acquisitions_info());

       if(!av.sorted())
           av.sort();

       auto sort_idx = av.get_kspace_order();
       auto sptr_enc = std::make_shared< RPEFourierEncoding > (RPEFourierEncoding());

       for(int i=0; i<sort_idx.size(); ++i)
       {
           sirf::AcquisitionsVector subset;
           av.get_subset(subset, sort_idx[i]);

           CFImage img;
           sptr_enc->backward(&img, subset);

           void* vptr_img = new CFImage(img);// god help me I don't trust this!

           ImageWrap iw(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, vptr_img);
           img_vec.append(iw);

           sptr_enc->forward(subset, &img);

           sptr_enc->backward(&img, subset);

           void* vptr_img_bfb = new CFImage(img);// god help me I don't trust this!
           ImageWrap iw_bfb(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, vptr_img_bfb);
           img_vec.append(iw_bfb);

           av.set_subset(subset, sort_idx[i]); //assume forward does not reorder the acquisitions
       }

       for(int i=0; i<img_vec.items(); ++i)
       {
           std::stringstream fname_output_img;
           fname_output_img << "output_" << __FUNCTION__ << "_image_" << i;
           write_cfimage_to_raw(fname_output_img.str(), img_vec.image_wrap(i));
       }

       std::stringstream fname_output_raw;
       fname_output_raw << "output_" << __FUNCTION__ << "_rawdata.h5";
       av.write(fname_output_raw.str());

       return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}




bool test_mracquisition_model_rpe_bwd(MRAcquisitionData& av)
{
    try
        {
            std::cout << "Running test " << __FUNCTION__ << std::endl;

            sirf::GRPETrajectoryPrep rpe_tp;
            rpe_tp.set_trajectory(av);

            sirf::GadgetronImagesVector img_vec;
            sirf::MRAcquisitionModel acquis_model;

            sirf::CoilSensitivitiesVector csm;
            csm.set_csm_smoothness(50);
            csm.calculate(av);

            auto sptr_encoder = std::make_shared<sirf::RPEFourierEncoding>(sirf::RPEFourierEncoding());
            acquis_model.set_encoder(sptr_encoder);

            acquis_model.bwd(img_vec, csm, av);

            for(int i=0; i<img_vec.items(); ++i)
            {
                std::stringstream fname_output;
                fname_output << "output_" << __FUNCTION__ << "_image_" << i;
                write_cfimage_to_raw(fname_output.str(), img_vec.image_wrap(i));
            }

            return true;

        }
        catch( std::runtime_error const &e)
        {
            std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
            std::cout << e.what() << std::endl;
            throw;
        }
}




int main ( int argc, char* argv[])
{

	try{

        std::string SIRF_PATH;
        if (argc==1)
            SIRF_PATH = getenv("SIRF_PATH");
        else
            SIRF_PATH = argv[1];

//        std::string data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian_Grappa2.h5";
        std::string data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian.h5";

        gadgetron::shared_ptr<MRAcquisitionData> sptr_ad(new AcquisitionsVector);
        AcquisitionsVector::set_as_template();

        AcquisitionsVector& av = (AcquisitionsVector&)*sptr_ad;

        av.read(data_path);

        sirf::preprocess_acquisition_data(av);
        av.sort();

        test_get_kspace_order(av);
        test_get_subset(av);

        test_GRPETrajectoryPrep_set_trajectory(av);

        test_CoilSensitivitiesVector_calculate(av);
        test_CoilSensitivitiesVector_get_csm_as_cfimage(av);

        test_bwd(av);

        std::string rpe_data_path = SIRF_PATH + "/data/examples/MR/Lowres_RPE_Dixon.h5";

        sirf::AcquisitionsVector rpe_av;
        rpe_av.read(rpe_data_path);

        sirf::preprocess_acquisition_data(rpe_av);
        rpe_av.sort();
        sirf::set_unit_dcf(rpe_av);


        test_get_rpe_trajectory(rpe_av);
        test_rpe_bwd(rpe_av);
        test_rpe_fwd(rpe_av);

        test_rpe_csm(rpe_av);

        test_mracquisition_model_rpe_bwd(rpe_av);

        test_number_of_encodings(sptr_ad);
        test_acq_mod_norm(sptr_ad);

        return 0;
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}





