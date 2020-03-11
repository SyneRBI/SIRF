#include <iostream>
#include <cstdlib>

#include "mrtest_auxiliary_funs.h"

#include "sirf/Gadgetron/encoding.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"
#include "sirf/Gadgetron/gadget_lib.h"

using namespace sirf;

bool test_get_kspace_order(void)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        std::string const fpath_input = "/media/sf_CCPPETMR/TestData/Input/xGadgetron/cGadgetron/";
        std::string fname_input = fpath_input + "CV_SR_64Cube_1Echo_10Dyn.h5";

        sirf::AcquisitionsVector av;
        av.read(fname_input);
        av.sort();

        auto kspace_sorting = av.get_kspace_order();

        fname_input = fpath_input + "CV_SR_128Cube_1Echo_3Dyn.h5";

        sirf::AcquisitionsVector av_contrast;
        av_contrast.read(fname_input);
        av_contrast.sort();

        auto kspace_sorting_contrast = av_contrast.get_kspace_order();


        fname_input = fpath_input + "CV_2D_Stack_144.h5";

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
        throw e;
    }
}

bool test_get_subset()
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        std::string const fpath_input = "/media/sf_CCPPETMR/TestData/Input/xGadgetron/cGadgetron/";
        std::string fname_input = fpath_input + "CV_SR_64Cube_1Echo_10Dyn.h5";

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
        throw e;
    }
}

int main ()
{
	try{

       test_get_kspace_order();
       test_get_subset();
        return 0;
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

