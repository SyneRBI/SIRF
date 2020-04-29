#include <iostream>
#include <cstdlib>

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"

#include "mrtest_auxiliary_funs.h"

using namespace sirf;

bool test_get_kspace_order(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av_slice;
        av_slice.read(fname_input);
        av_slice.sort();

        auto kspace_sorting_slice = av_slice.get_kspace_order();

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_get_subset(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);

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

bool test_CoilSensitivitiesVector_calculate(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        sirf::preprocess_acquisition_data(av);

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

bool test_CoilSensitivitiesVector_get_csm_as_cfimage(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        sirf::preprocess_acquisition_data(av);

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



int main ( int argc, char* argv[])
{

	try{

        std::string SIRF_PATH;
        if (argc==1)
            SIRF_PATH = getenv("SIRF_PATH");
        else
            SIRF_PATH = argv[1];

        std::string data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian_Grappa2.h5";

//        test_get_kspace_order(data_path);
//        test_get_subset(data_path);
        test_CoilSensitivitiesVector_calculate(data_path);
        test_CoilSensitivitiesVector_get_csm_as_cfimage(data_path);
        return 0;
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

