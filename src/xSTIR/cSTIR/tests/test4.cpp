#include <iostream>
#include <cstdlib>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"
USING_NAMESPACE_STIR
USING_NAMESPACE_ECAT

#include "cstir.h"
#include "object.h"
#include "stir_x.h"
//#include "SIRF/common/envar.h"

int main()
{
    try{

			//std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
			std::string SIRF_path = std::getenv("SIRF_PATH");
			if (SIRF_path.length() < 1) {
            std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
            return 1;
        }

        std::string path = SIRF_path + "/data/examples/PET/mMR/";

        string f_listmode   = path + "list.l.hdr";
        string f_template   = path + "mMR_template_span11_small.hs";

        // Listmode to sinograms
        ListmodeToSinograms converter;
        converter.set_input   ( f_listmode  );
        converter.set_output  ( "proj_data" );
        converter.set_template( f_template  );
        converter.set_time_interval(0,10);
        converter.set_up();
        converter.estimate_randoms();
		return 0;
	}
    catch (...)
    {
		return 1;
	}
}