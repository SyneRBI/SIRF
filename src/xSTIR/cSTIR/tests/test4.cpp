#include <iostream>
#include <cstdlib>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"

#include "cstir.h"
#include "object.h"

#include "sirf/cSTIR/stir_x.h"
#include "sirf/common/getenv.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

int test4()
{
	try {

		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

		std::string path = SIRF_path + "/data/examples/PET/mMR/";

		string f_listmode = path + "list.l.hdr";
		string f_template = path + "mMR_template_span11_small.hs";

		// Listmode to sinograms
		ListmodeToSinograms converter;
		converter.set_input(f_listmode);
		converter.set_output("proj_data");
		converter.set_template(f_template);
		converter.set_time_interval(0, 10);
		converter.set_up();
		converter.estimate_randoms();
		return 0;
	}
	catch (...)
	{
		return 1;
	}
}

//int test5();

int main()
{
	return test4();
	//return test5();
}
