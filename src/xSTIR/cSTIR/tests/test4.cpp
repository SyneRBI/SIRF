#include <iostream>

#include "cstir.h"
#include "object.h"
#include "stir_types.h"
#include "stir_data_containers.h"
#include "stir_x.h"
#include "SIRF/common/envar.h"

using namespace std;

int main()
{
    try{

        std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
        if (SIRF_path.length() < 1) {
            std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
            return 1;
        }

        std::string path = SIRF_path + "/data/examples/PET/";

        string f_listmode   = path + "list.l.hdr.STIR";
        string f_template   = path + "template_span11.hs";

        // Listmode to sinograms
        ListmodeToSinograms converter;
        converter.set_input   ( f_listmode  );
        converter.set_output  ( "proj_data" );
        converter.set_template( f_template  );
        converter.set_time_interval(0,10);
        converter.set_up();
        converter.estimate_randoms();
    }
    catch (...)
    {
    }
    return 0;
}