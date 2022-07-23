#include "tests_c_interface.h"

#include "test_input_filenames.h"
#include "auxiliary_testing_functions.h"

#include "sirf/iUtilities/DataHandle.h"
#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/cDynamicSimulation/cdynamicsimulation.h"

using namespace sirf;

bool test_simulation_interface::test_bin_data_from_handle( void )
{

    std::cout << "--- Running "<< __FUNCTION__ << std::endl;

	try
	{
		bool test_succesful = true;

		AcquisitionsVector acq_vec;
		std::string fpath_testdata("/media/sf_CCPPETMR/TestData/Input/xDynamicSimulation/cDynamicSimulation/TemplateData/MR/CV_nav_cart_64Cube_1Echo.h5");

		acq_vec.read(fpath_testdata);

        std::shared_ptr<MRAcquisitionData> 
			sptr_ad(new AcquisitionsVector(acq_vec));

		void* ptr_ad =  newObjectHandle<MRAcquisitionData>(sptr_ad);

		int const num_bins = 4;
		MRMotionDynamic motion_dyn(num_bins);
		SignalContainer mock_signal = aux_test::get_generic_cardiac_signal(acq_vec);
		motion_dyn.set_dynamic_signal(mock_signal);

        std::shared_ptr<MRMotionDynamic> 
			sptr_motiondyn(new MRMotionDynamic(motion_dyn));

		void* ptr_dyn =  newObjectHandle<MRMotionDynamic>(sptr_motiondyn);
        
		cDS_setMRAcquisitions(ptr_dyn, ptr_ad);

		return true;
    }
	catch( std::runtime_error const &e)
	{
		cout << "Exception caught " <<__FUNCTION__ <<" .!" <<endl;
		cout << e.what() << endl;
		throw e;
	}
}