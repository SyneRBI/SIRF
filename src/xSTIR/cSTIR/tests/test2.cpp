#include <string>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"
USING_NAMESPACE_STIR
USING_NAMESPACE_ECAT

#include "cstir.h"
#include "handle.h"
#include "stir_types.h"
//#include "SIRF/common/envar.h"

void* TMP_HANDLE;

#define GET_FLOAT(V, F) \
	TMP_HANDLE = F; \
	if (execution_status(TMP_HANDLE)) break; \
	V = floatDataFromHandle(TMP_HANDLE); \
	deleteDataHandle(TMP_HANDLE)

int test2()
{
	std::string filename;
	int dim[3];
	float at_value = 0.05f*0;
	float bt_value = 0.1f*0;
	float s, t;
	void* handle = 0;
	void* matrix = 0;
	void* image = 0;
	void* img = 0;
	void* am = 0;
	void* ad = 0;
	void* fd = 0;
	void* at = 0;
	void* bt = 0;
	void* nd = 0;
	void* norm = 0;
	void* prior = 0;
	void* obj_fun = 0;
	void* filter = 0;
	void* recon = 0;
	void* diff = 0;
	void* sm = 0;
	void* ai = 0;

	//std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
	std::string SIRF_path = std::getenv("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
		return 1;
	}
	std::string path = SIRF_path + "/data/examples/PET/";

	TextWriter w;
	openChannel(0, &w);

	for (;;) {
		//filename = SIRF_path + "/examples/Python/PET/my_image.hv";
		//HANDLE(image, cSTIR_objectFromFile("Image", filename.c_str()));
		//cSTIR_getImageDimensions(image, (size_t)&dim[0]);
		//std::cout << "image dimensions: " 
		//	<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

		HANDLE(matrix, cSTIR_newObject("RayTracingMatrix"));
		CALL(cSTIR_setParameter
			(matrix, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2)));

		//filename = path + "my_forward_projection.hs";
		filename = "sinograms_f1g1d0b0.hs";
		HANDLE(ad, cSTIR_objectFromFile("AcquisitionData", filename.c_str()));
		cSTIR_getAcquisitionsDimensions(ad, (size_t)&dim[0]);
		std::cout << "acquisition data dimensions: "
			<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

		HANDLE(image, cSTIR_imageFromAcquisitionData(ad));
		cSTIR_getImageDimensions(image, (size_t)&dim[0]);
		std::cout << "image dimensions: " 
			<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

		//HANDLE(at, cSTIR_acquisitionsDataFromTemplate(ad));
		//cSTIR_getAcquisitionsDimensions(at, (size_t)&dim[0]);
		//cSTIR_fillAcquisitionsData(at, at_value);
		//HANDLE(bt, cSTIR_acquisitionsDataFromTemplate(ad));
		//cSTIR_getAcquisitionsDimensions(bt, (size_t)&dim[0]);
		//cSTIR_fillAcquisitionsData(bt, bt_value);
		//HANDLE(nd, cSTIR_acquisitionsDataFromTemplate(ad));
		//cSTIR_getAcquisitionsDimensions(nd, (size_t)&dim[0]);
		//cSTIR_fillAcquisitionsData(nd, 2.0);

		HANDLE(am, cSTIR_newObject("AcqModUsingMatrix"));
		//CALL(cSTIR_setParameter(am, "AcquisitionModel", "additive_term", at));
		//CALL(cSTIR_setParameter(am, "AcquisitionModel", "normalisation", nd));
		CALL(cSTIR_setParameter(am, "AcqModUsingMatrix", "matrix", matrix));
		CALL(cSTIR_setupAcquisitionModel(am, ad, image));

		filename = path + "mu_map.hv";
		HANDLE(ai, cSTIR_objectFromFile("Image", filename.c_str()));
		HANDLE(sm, cSTIR_createPETAttenuationModel(ai, am));
		CALL(cSTIR_setParameter(am, "AcquisitionModel", "asm", sm));
		CALL(cSTIR_setupAcquisitionSensitivityModel(sm, ad));
		std::cout << "ok\n";

		//std::cout << "projecting...\n";
		//HANDLE(fd, cSTIR_acquisitionModelFwd(am, image));
		//cSTIR_getAcquisitionsDimensions(fd, (size_t)&dim[0]);
		//std::cout << "simulated acquisition data dimensions: "
		//	<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		//GET_FLOAT(s, cSTIR_norm(ad));
		//GET_FLOAT(t, cSTIR_norm(fd));
		//diff = cSTIR_axpby(1/s, ad, -1/t, fd);
		//GET_FLOAT(s, cSTIR_norm(diff));
		//deleteDataHandle(diff);
		//std::cout << "acq diff: " << s << '\n';

		//HANDLE(img, cSTIR_acquisitionModelBwd(am, fd));
		//cSTIR_getImageDimensions(img, (size_t)&dim[0]);
		//std::cout << "backprojected image dimensions: " 
		//	<< dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		//if (at_value == 0 && bt_value == 0) {
		//	GET_FLOAT(s, cSTIR_dot(img, image));
		//	std::cout << s << " = " << t*t << '\n';
		//}
		//deleteDataHandle(img);

		HANDLE(prior, cSTIR_newObject("QuadraticPrior"));
		
		HANDLE(prior2, cSTIR_newObject("PLSPrior"));
		handle = floatDataHandle(0.5);
		CALL(cSTIR_setParameter(prior2, "PLSPrior", "penalisation_factor", handle));
		deleteDataHandle(handle);
		handle = floatDataHandle(0.5);
		CALL(cSTIR_setParameter(prior2, "PLSPrior", "alpha", handle));
		deleteDataHandle(handle);
		handle = floatDataHandle(0.5);
		CALL(cSTIR_setParameter(prior2, "PLSPrior", "eta", handle));
		deleteDataHandle(handle);
		handle = charDataHandle("test");
		CALL(cSTIR_setParameter(prior2, "PLSPrior", "kappa_filename", handle));
		deleteDataHandle(handle);
		handle = charDataHandle("test");
		CALL(cSTIR_setParameter(prior2, "PLSPrior", "anatomical_filename", handle));
		deleteDataHandle(handle);
		cSTIR_setupPLSPrior(prior2);

		std::string obj_fun_name
			("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
		HANDLE(obj_fun, cSTIR_newObject(obj_fun_name.c_str()));
		CALL(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "acquisition_model", am));
		//CALL(cSTIR_setParameter
		//	(obj_fun, obj_fun_name.c_str(), "acquisition_data", fd));
		CALL(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "acquisition_data", ad));
		handle = charDataHandle("true");
		CALL(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "zero_seg0_end_planes", handle));
		deleteDataHandle(handle);
		int max_seg_num = 4; // causes crash if < 4
		handle = intDataHandle(max_seg_num);
		CALL(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "max_segment_num_to_process", handle));
		deleteDataHandle(handle);
		CALL(cSTIR_setParameter
			(obj_fun, "GeneralisedObjectiveFunction", "prior", prior));
		//CALL(cSTIR_setupObjectiveFunction(obj_fun, image));

		std::cout << "ok\n";

		HANDLE(filter, cSTIR_newObject("TruncateToCylindricalFOVImageProcessor"));

		int num_subiterations = 2;
		HANDLE(recon, cSTIR_objectFromFile("OSMAPOSLReconstruction", ""));
		handle = charDataHandle("reconstructedImage");
		CALL(cSTIR_setParameter
			(recon, "Reconstruction", "output_filename_prefix", handle));
		deleteDataHandle(handle);
		handle = intDataHandle(12);
		CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "num_subsets", handle));
		deleteDataHandle(handle);
		handle = intDataHandle(num_subiterations);
		CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "num_subiterations", handle));
		CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "save_interval", handle));
		deleteDataHandle(handle);
		handle = intDataHandle(1);
		CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "inter_iteration_filter_interval", 
			handle));
		deleteDataHandle(handle);
		CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "objective_function", obj_fun));
		CALL(cSTIR_setParameter
			(recon, "IterativeReconstruction", "inter_iteration_filter_type", filter));
		handle = charDataHandle("multiplicative");
		CALL(cSTIR_setParameter(recon, "OSMAPOSL", "MAP_model", handle));
		deleteDataHandle(handle);
		std::cout << "ok\n";
		CALL(cSTIR_setupReconstruction(recon, image));

		for (int iter = 0; iter < num_subiterations; iter++) {
			std::cout << "iteration " << iter << '\n';
			cSTIR_updateReconstruction(recon, image);
		}

		break;
	}

	deleteDataHandle(image);
	deleteDataHandle(img);
	deleteDataHandle(ad);
	deleteDataHandle(am);
	deleteDataHandle(matrix);
	deleteDataHandle(at);
	deleteDataHandle(bt);
	deleteDataHandle(nd);
	deleteDataHandle(fd);
	deleteDataHandle(recon);
	deleteDataHandle(filter);
	deleteDataHandle(prior);
	deleteDataHandle(obj_fun);

	return 0;
}

