#include <fstream>
#include <string>

//#include "stir/listmode/LmToProjData.h"

#include "cstir.h"
#include "object.h"
#include "stir_types.h"
#include "stir_data_containers.h"
#include "stir_x.h"
#include "SIRF/common/envar.h"

int test3()
{
	std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
		return 1;
	}
	std::string path = SIRF_path + "/data/examples/PET/";
	std::string filename = path + "list.l.hdr.STIR";

	//LmToProjData lm_data("lm_to_projdata.par");
	//lm_data.process_data();
	ListmodeToSinograms converter; // ("lm_to_projdata.par");
	//converter.set_input("list.l.hdr.STIR");
	converter.set_input(filename);
	converter.set_output("proj_data");
	filename = path + "template_span11.hs";
	converter.set_template(filename);
	//converter.set_template("template_span11.hs");
	converter.set_time_interval(0, 10);
	converter.set_up();
	//converter.process_data();
	converter.compute_fan_sums();
	converter.compute_singles();
	converter.estimate_randoms();

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
		handle = cSTIR_convertListmodeToSinograms(lm2s);
		deleteDataHandle(handle);
		//CALL(cSTIR_convertListmodeToSinograms(lm2s));

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
		void* sm = 0;
		for (;;) {
			void* h = charDataHandle(input_filename.c_str());
			HANDLE(sm, cSTIR_createPETAcquisitionSensitivityModel(h, "n"));
			break;
		}
		deleteDataHandle(sm);
		//shared_ptr<PETAcquisitionSensitivityModel> 
		//	sptr(new PETAcquisitionSensitivityModel(input_filename));
		//shared_ptr<BinNormalisationFromECAT8>
		//	sptr_n(new BinNormalisationFromECAT8(input_filename));
		//shared_ptr<CListModeData> lm_data_ptr;
		//lm_data_ptr = stir::read_from_file<CListModeData>(input_filename);
	}
	catch (...)
	{
	}
	return 0;
}

int test3c()
{
	try {
		std::string filename ="mu_map.hv";
		void* image = 0;
		void* sm = 0;
		for (;;) {
			HANDLE(image, cSTIR_objectFromFile("Image", filename.c_str()));
			HANDLE(sm, cSTIR_createPETAcquisitionSensitivityModel(image, "i"));
			break;
		}
		deleteDataHandle(image);
		deleteDataHandle(sm);
		//shared_ptr<PETImageData> sptr_img(new PETImageData(filename));
		//shared_ptr<PETAcquisitionSensitivityModel>
		//	sptr(new PETAcquisitionSensitivityModel(*sptr_img));
		//shared_ptr<BinNormalisationFromAttenuationImage>
		//	sptr_n(new BinNormalisationFromAttenuationImage(sptr_img->data_sptr()));
	}
	catch (...)
	{
	}
	return 0;
}

int test3d()
{
	try {
		void* ad = 0;
		void* sm = 0;
		for (;;) {
			int span = 1;
			int max_ring_diff = -1;
			int view_mash_factor = 1;
			std::cout << "creating acquisition data...";
			HANDLE(ad, cSTIR_acquisitionsDataFromScannerInfo
				("Siemens_mMR", span, max_ring_diff, view_mash_factor));
			std::cout << "ok\n";
			CALL(cSTIR_fillAcquisitionsData(ad, 1.0f));
			std::cout << "creating acquisition sensitivity data...";
			HANDLE(sm, cSTIR_createPETAcquisitionSensitivityModel(ad, "s"));
			std::cout << "ok\n";
			break;
		}
		deleteDataHandle(ad);
		deleteDataHandle(sm);
	}
	catch (...)
	{
	}
	return 0;
}