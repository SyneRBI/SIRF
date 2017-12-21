#include <fstream>
#include <string>

//#include "stir/listmode/LmToProjData.h"

#include "object.h"
#include "stir_types.h"
#include "stir_data_containers.h"
#include "stir_x.h"
#include "envar.h"

int test3()
{
	//LmToProjData lm_data("lm_to_projdata.par");
	//lm_data.process_data();
	ListmodeToSinograms converter; // ("lm_to_projdata.par");
	converter.set_input("list.l.hdr.STIR");
	converter.set_output("proj_data");
	converter.set_template("template_span11.hs");
	converter.set_time_interval(0, 10);
	converter.set_up();
	converter.process_data();

	return 0;
}