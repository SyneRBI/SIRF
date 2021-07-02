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
\ingroup MR
\brief MR related C++ tests

\author Johannes Mayer
\author SyneRBI
*/

#include <iostream>
#include <cstdlib>
#include <numeric>
#include <vector>
#include <random>

#include <ismrmrd/xml.h>

#include "sirf/Gadgetron/chain_lib.h"

#include "sirf/common/DataContainer.h"
#include "sirf/common/getenv.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"

#include "sirf/Gadgetron/FourierEncoding.h"
#include "sirf/Gadgetron/TrajectoryPreparation.h"

#include "mrtest_auxiliary_funs.h"

using namespace sirf;

bool const mr_cpp_tests_writefiles = true;

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

bool test_ISMRMRDImageData_from_MRAcquisitionData(MRAcquisitionData& av)
{
     try
    {
        using namespace ISMRMRD;

        std::cout << "Running test " << __FUNCTION__ << std::endl;

        GadgetronImagesVector iv(av);
        
        bool test_successful = true;

        //
        int const num_images = iv.number();
        int const num_kspace_dims = av.get_kspace_sorting().size();
        test_successful *= (num_images == num_kspace_dims);


        //
        std::vector<Encoding> enc_vec = av.acquisitions_info().get_IsmrmrdHeader().encoding;
        Encoding enc = enc_vec[0];
        EncodingSpace rec_space = enc.reconSpace;

        MatrixSize rawdata_recon_matrix = rec_space.matrixSize;
        FieldOfView_mm rawdata_recon_FOV = rec_space.fieldOfView_mm;

        for(int i=0; i<iv.number(); ++i)
        {
            CFImage* const ptr_img = (CFImage*)(iv.sptr_image_wrap(i)->ptr_image());            

            test_successful *= (ptr_img->getMatrixSizeX() == rawdata_recon_matrix.x);
            test_successful *= (ptr_img->getMatrixSizeY() == rawdata_recon_matrix.y);
            test_successful *= (ptr_img->getMatrixSizeZ() == rawdata_recon_matrix.z);

        	test_successful *= (ptr_img->getFieldOfViewX() == rawdata_recon_FOV.x);
            test_successful *= (ptr_img->getFieldOfViewY() == rawdata_recon_FOV.y);
            test_successful *= (ptr_img->getFieldOfViewZ() == rawdata_recon_FOV.z);
        }

        if(test_successful)
            return test_successful;
        else{
            throw std::runtime_error("The test for images from acquisition data failed.");                
        }

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_ISMRMRDImageData_reorienting(MRAcquisitionData& av)
{
    try
    {
        using namespace ISMRMRD;

        std::cout << "Running test " << __FUNCTION__ << std::endl;

        GadgetronImagesVector iv(av);
        std::shared_ptr<const VoxelisedGeometricalInfo3D > sptr_geom_info = iv.get_geom_info_sptr();

        VoxelisedGeometricalInfo3D::Size input_size = sptr_geom_info->get_size();

        std::default_random_engine generator;
        generator.seed(100000); // 
        float const maximum_random_range = 100.f;
        std::uniform_real_distribution<float> distribution(0, maximum_random_range);

        float const rand_angle = distribution(generator);
        std::cout << "Rotating by an angle of : " << rand_angle << std::endl;
        VoxelisedGeometricalInfo3D::DirectionMatrix random_rotation{
            std::cos(rand_angle), -std::sin(rand_angle), 0,
            std::sin(rand_angle),  std::cos(rand_angle), 0,
            0, 0, 1};

        VoxelisedGeometricalInfo3D::DirectionMatrix direction = sptr_geom_info->get_direction();
        VoxelisedGeometricalInfo3D::DirectionMatrix rotated_direction = direction;
        
        for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        {   
            rotated_direction[i][j] = 0.0;

            for(int k=0; k<3; ++k)
            for(int l=0; l<3; ++l)
                rotated_direction[i][j] += random_rotation[i][k]*direction[k][l]*random_rotation[j][l];

        }
        VoxelisedGeometricalInfo3D::Offset offset{distribution(generator),distribution(generator),distribution(generator)};
        VoxelisedGeometricalInfo3D::Spacing spacing{distribution(generator),distribution(generator),distribution(generator)};                    

        VoxelisedGeometricalInfo3D random_new_geometry(offset, spacing, input_size, rotated_direction);

        if(!iv.can_reorient(*sptr_geom_info, random_new_geometry, false))
            throw std::runtime_error("You can not perform this arbitrary reorienting.");

        iv.reorient(random_new_geometry);
        
        bool test_successful = true;

        if(test_successful)
            return test_successful;
        else{
            throw std::runtime_error("The test for images from acquisition data failed.");                
        }

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

        if(mr_cpp_tests_writefiles)
            sirf::write_imagevector_to_raw(__FUNCTION__, csv);

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

            if(mr_cpp_tests_writefiles)
            {
                std::stringstream fname_out;
                fname_out << "output_" << __FUNCTION__ << "_" << i;
                sirf::write_cfimage_to_raw(fname_out.str(), img);
            }
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

bool test_acq_mod_adjointness(MRAcquisitionData& ad)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        std::cout << "Assessing if operator E is adjoint by comparing <Eh k, i> and <k, E i> <>" << std::endl;        
        std::cout << "where E is the MR acquisition model, k is k-space data and i is a random complex image." << std::endl;        

        // setup the acquisition model                 
        sirf::MRAcquisitionModel AM = sirf::get_prepared_MRAcquisitionModel(ad);
        
        auto sptr_bwd = AM.bwd(ad);

        sirf::Dimensions dims = sptr_bwd->dimensions();

        // generate a random image to project onto
        int const num_total_pixels = dims["x"]*dims["y"]*dims["z"]*dims["c"]*dims["n"];

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(0.0,1.0);

        std::vector<complex_float_t> random_data;
        for(int i=0; i<num_total_pixels; ++i)
        {
            float const real_part = distribution(generator);
            float const imag_part = distribution(generator);
            complex_float_t curr_number = std::complex<float>(real_part, imag_part);
            random_data.push_back(curr_number);
        }

        std::shared_ptr<ISMRMRDImageData> sptr_random = std::move(sptr_bwd->clone());
        sptr_random->set_data(&random_data[0]);

        complex_float_t Eh_kdat_Dot_img;
        sptr_bwd->dot(*sptr_random, &Eh_kdat_Dot_img);

        auto sptr_fwd = AM.fwd(*sptr_random);

        complex_float_t E_img_Dot_kdat;
        ad.dot(*sptr_fwd, &E_img_Dot_kdat);

        std::cout << "Backward kdata dot random image: " << Eh_kdat_Dot_img << std::endl;
        std::cout << "Forward random image dot kdata : " << E_img_Dot_kdat  << std::endl;

        float const order_of_magnitude = std::min( std::abs(Eh_kdat_Dot_img), std::abs(E_img_Dot_kdat));
        float const diff_in_scalar_prod = std::abs(Eh_kdat_Dot_img - E_img_Dot_kdat);
        float const tolerance = 0.0001;
        
        std::cout <<"Level of non-adjointness is given by: |<Eh k, i> - <k, E i> <>|/max(|<Eh k, i>,<k, E i> <>|) = " <<  diff_in_scalar_prod/order_of_magnitude << std::endl;
        std::cout <<"Accepting a relative tolerance of " << tolerance << std::endl;
        
        bool ok = diff_in_scalar_prod/order_of_magnitude < tolerance;
    
        return ok;
    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_acq_mod_norm(gadgetron::shared_ptr<MRAcquisitionData> sptr_ad)
{

try
    {   
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        MRAcquisitionModel AM = sirf::get_prepared_MRAcquisitionModel(*sptr_ad);

        std::cout << "Back projection ... "<< std::endl;
        auto sptr_bwd = AM.bwd(*sptr_ad);
        std::cout << "Forward projection ... "<< std::endl;
        auto sptr_fwd = AM.fwd(*sptr_bwd);

        float im_norm = sptr_bwd->norm();
        float sd_norm = sptr_fwd->norm();

        AM.set_up(sptr_fwd, sptr_bwd);
        
        std::cout << "AM norm calculation ... " << std::endl;
        float am_norm = AM.norm();

        std::cout << "\nChecking the acquisition model norm:\n";
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
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
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

        if(mr_cpp_tests_writefiles)
            sirf::write_imagevector_to_raw(__FUNCTION__, img_vec);

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}


bool test_TrajectoryPreparation_constructors( void )
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

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



bool test_set_rpe_trajectory(AcquisitionsVector av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;
        sirf::GRPETrajectoryPrep rpe_tp;

        rpe_tp.set_trajectory(av);

        if(mr_cpp_tests_writefiles)
        {
            std::stringstream fname_output;
            fname_output << "output_" << __FUNCTION__ << ".h5";
            av.write(fname_output.str());
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

#ifdef GADGETRON_TOOLBOXES_AVAILABLE
#warning "INCLUDING THE RADIAL TESTS INTO THE C++ TESTS"
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

       if(mr_cpp_tests_writefiles)
           sirf::write_imagevector_to_raw(__FUNCTION__, csm);
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

           CFImage* ptr_img = new CFImage();// god help me I don't trust this!
           ImageWrap iw(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, ptr_img);

           sptr_enc->backward(*ptr_img, subset);

           

           img_vec.append(iw);

       }

       if(mr_cpp_tests_writefiles)
           sirf::write_imagevector_to_raw(__FUNCTION__, img_vec);


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

           CFImage* ptr_img = new CFImage();
           ImageWrap iw(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, ptr_img);

           sptr_enc->backward(*ptr_img, subset);

           img_vec.append(iw);

           sptr_enc->forward(subset, *ptr_img);

           sptr_enc->backward(*ptr_img, subset);

           CFImage* ptr_img_bfb = new CFImage(*ptr_img);// god help me I don't trust this!
           ImageWrap iw_bfb(ISMRMRD::ISMRMRD_DataTypes::ISMRMRD_CXFLOAT, ptr_img_bfb);
           
           img_vec.append(iw_bfb);

           av.set_subset(subset, sort_idx[i]); //assume forward does not reorder the acquisitions
       }

       if(mr_cpp_tests_writefiles)
       {
           sirf::write_imagevector_to_raw(__FUNCTION__, img_vec);

           std::stringstream fname_output_raw;
           fname_output_raw << "output_" << __FUNCTION__ << "_rawdata.h5";
           av.write(fname_output_raw.str());
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

            if(mr_cpp_tests_writefiles)
                sirf::write_imagevector_to_raw(__FUNCTION__, img_vec);

            return true;

        }
        catch( std::runtime_error const &e)
        {
            std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
            std::cout << e.what() << std::endl;
            throw;
        }
}
#endif


int main ( int argc, char* argv[])
{
	try{
		
        std::string SIRF_PATH;
        if (argc==1)
            SIRF_PATH = sirf::getenv("SIRF_PATH", true);
        else
            SIRF_PATH = argv[1];

//        std::string data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian_Grappa2.h5";
        std::string data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian.h5";

        std::shared_ptr<MRAcquisitionData> sptr_ad(new AcquisitionsVector);
        AcquisitionsVector& av = (AcquisitionsVector&)*sptr_ad;
        av.read(data_path);

        sirf::preprocess_acquisition_data(av);
        av.sort();

        bool ok = true;

        ok *= test_get_kspace_order(av);
        ok *= test_get_subset(av);

        ok *= test_ISMRMRDImageData_from_MRAcquisitionData(av);
        ok *= test_ISMRMRDImageData_reorienting(av);

        ok *= test_CoilSensitivitiesVector_calculate(av);
        ok *= test_CoilSensitivitiesVector_get_csm_as_cfimage(av);

        ok *= test_bwd(av);

        ok *= test_acq_mod_adjointness(av);
        ok *= test_acq_mod_norm(sptr_ad);


        #ifdef GADGETRON_TOOLBOXES_AVAILABLE
        #warning "RUNNING THE RADIAL TESTS FOR C++."
            std::string rpe_data_path = SIRF_PATH + "/data/examples/MR/zenodo/3D_RPE_Lowres.h5";
            sirf::AcquisitionsVector rpe_av;
            rpe_av.read(rpe_data_path);



            sirf::preprocess_acquisition_data(rpe_av);
            rpe_av.sort();
            sirf::set_unit_dcf(rpe_av);


            ok *= test_set_rpe_trajectory(rpe_av);
            ok *= test_rpe_bwd(rpe_av);
            ok *= test_rpe_fwd(rpe_av);

            ok *= test_rpe_csm(rpe_av);

            ok *= test_mracquisition_model_rpe_bwd(rpe_av);
            ok *= test_acq_mod_adjointness(rpe_av);

            auto sptr_rpe_av = std::make_shared<AcquisitionsVector>(rpe_av);
            sirf::GRPETrajectoryPrep rpe_tp;
            rpe_tp.set_trajectory(*sptr_rpe_av);
            ok *= test_acq_mod_norm(sptr_rpe_av);
        #endif



        if(ok)
            return 0;
        else
        { 
            std::cerr << "The code ran but some quantitative tests must have failed"<<std::endl;
            return EXIT_FAILURE;
        }
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}





