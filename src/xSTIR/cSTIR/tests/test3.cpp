#include <fstream>
#include <string>

//#include "stir/listmode/LmToProjData.h"

#include "object.h"
#include "stir_types.h"
#include "stir_data_containers.h"
#include "stir_x.h"
#include "SIRF/common/envar.h"

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

#include "cstir.h"
#include "handle.h"

int test3a()
{
	std::string filename;
	void* lm2s = 0;
	//void* h = 0;
	void* handle = 0;
	float interval[] = { 0, 10 };

	std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
		return 1;
	}
	std::string path = SIRF_path + "/data/examples/PET/";
	filename = path + "list.l.hdr.STIR";

	for (;;) {
		//HANDLE(lm2s, cSTIR_objectFromFile("ListmodeToSinograms", "lm_to_projdata.par"));
		HANDLE(lm2s, cSTIR_newObject("ListmodeToSinograms"));
		//handle = charDataHandle("list.l.hdr.STIR");
		handle = charDataHandle(filename.c_str());
		CALL(cSTIR_setParameter
			(lm2s, "ListmodeToSinograms", "input", handle));
		deleteDataHandle(handle);
		handle = charDataHandle("proj_data");
		CALL(cSTIR_setParameter
			(lm2s, "ListmodeToSinograms", "output", handle));
		deleteDataHandle(handle);
		filename = path + "template_span11.hs";
		handle = charDataHandle(filename.c_str());
		//handle = charDataHandle("template_span11.hs");
		CALL(cSTIR_setParameter
			(lm2s, "ListmodeToSinograms", "template", handle));
		deleteDataHandle(handle);
		CALL(cSTIR_setListmodeToSinogramsInterval(lm2s, (size_t)interval));
		CALL(cSTIR_setupListmodeToSinogramsConverter(lm2s));
		CALL(cSTIR_convertListmodeToSinograms(lm2s));

		break;
	}

	deleteDataHandle(lm2s);

	return 0;
}

int test3b()
{
	try {
		std::string input_filename =
			"C:/Users/wps46139/Documents/GitHub/SIRF/data/examples/PET/norm.n.hdr.STIR";
		//"C:/Users/wps46139/Documents/GitHub/SIRF/data/examples/PET/list.l.hdr.STIR";
		shared_ptr<BinNormalisationFromECAT8>
			sptr_n(new BinNormalisationFromECAT8(input_filename));
		//shared_ptr<CListModeData> lm_data_ptr;
		//lm_data_ptr = stir::read_from_file<CListModeData>(input_filename);
	}
	catch (...)
	{
	}
}