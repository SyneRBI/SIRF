#include <string>

//#include "stir/TextWriter.h"

#include "cstir.h"
#include "data_handle.h"
#include "iutilities.h"
#include "stir.h"
#include "tests.h"

void* TMP_HANDLE;

#define RUN(F) TMP_HANDLE = F; if (execution_status(TMP_HANDLE, 1)) break

int execution_status(void* handle, int clear = 0)
{
	int s = executionStatus(handle);
	if (s)
		std::cout << executionError(handle) << '\n';
	if (clear)
		deleteDataHandle(handle);
	return s;
}

void test2()
{
	std::string path("../../examples/");
	std::string filename;
	int status;
	int dim[3];
	void* h = 0;
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

	TextWriter w;
	openChannel(0, &w);

	for (;;) {
		image = cSTIR_objectFromFile("Image", (path + "my_image.hv").c_str());
		cSTIR_getImageDimensions(image, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		size_t image_size = dim[0] * dim[1] * dim[2];
		double* image_data = new double[image_size];
		cSTIR_getImageData(image, (size_t)image_data);

		matrix = cSTIR_newObject("RayTracingMatrix");
		handle = cSTIR_setParameter
			(matrix, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2));
		//status = execution_status(handle);
		//if (status)
		//	break;
		ad = cSTIR_objectFromFile
			("AcquisitionData", (path + "my_forward_projection.hs").c_str());
		status = execution_status(ad);
		if (status)
			break;
		cSTIR_getAcquisitionsDimensions(ad, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		size_t size = dim[0] * dim[1] * dim[2];

		at = cSTIR_acquisitionsDataFromTemplate(ad);
		cSTIR_getAcquisitionsDimensions(at, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		cSTIR_fillAcquisitionsData(at, 0.05);
		bt = cSTIR_acquisitionsDataFromTemplate(ad);
		cSTIR_getAcquisitionsDimensions(bt, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		cSTIR_fillAcquisitionsData(bt, 0.05);
		nd = cSTIR_acquisitionsDataFromTemplate(ad);
		cSTIR_getAcquisitionsDimensions(nd, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		cSTIR_fillAcquisitionsData(nd, 2.0);

		am = cSTIR_newObject("AcqModUsingMatrix");
		RUN(cSTIR_setParameter(am, "AcqModUsingMatrix", "matrix", matrix));
		RUN(cSTIR_setupAcquisitionModel(am, ad, image));
		fd = cSTIR_acquisitionModelFwd(am, image, "");
		cSTIR_getAcquisitionsDimensions(fd, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

		//double* adv = new double[size];
		//double* fdv = new double[size];

		//cSTIR_getAcquisitionsData(ad, (size_t)adv);
		//cSTIR_getAcquisitionsData(fd, (size_t)fdv);
		//std::cout << "acq diff: " << diff(size, adv, fdv) << '\n';

		//img = cSTIR_acquisitionModelBwd(am, fd);
		//cSTIR_getImageDimensions(img, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		//size_t im_size = dim[0] * dim[1] * dim[2];
		//double* img_data = new double[im_size];
		//cSTIR_getImageData(img, (size_t)img_data);

		//std::cout << dot(size, adv, fdv) << '\n';
		//std::cout << dot(im_size, image_data, img_data) << '\n';

		//delete[] adv;
		//delete[] fdv;
		//delete[] img_data;

		//RUN(cSTIR_setParameter(am, "AcquisitionModel", "additive_term", at));
		//RUN(cSTIR_setParameter(am, "AcquisitionModel", "normalisation", nd));
		//deleteDataHandle(fd);
		//fd = cSTIR_acquisitionModelFwd(am, image, "");

		prior = cSTIR_newObject("QuadraticPrior");

		std::string obj_fun_name
			("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
		obj_fun = cSTIR_newObject(obj_fun_name.c_str());
		RUN(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "acquisition_model", am));
		RUN(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "proj_data_sptr", fd));
		handle = charDataHandle("true");
		RUN(cSTIR_setParameter
			(obj_fun, obj_fun_name.c_str(), "zero_seg0_end_planes", handle));
		deleteDataHandle(handle);
		// causes crush
		//handle = intDataHandle(3);
		//RUN(cSTIR_setParameter
		//	(obj_fun, obj_fun_name.c_str(), "max_segment_num_to_process", handle));
		//deleteDataHandle(handle);
		RUN(cSTIR_setParameter
			(obj_fun, "GeneralisedObjectiveFunction", "prior", prior));

		filter = cSTIR_newObject("TruncateToCylindricalFOVImageProcessor");

		int num_subiterations = 2;
		recon = cSTIR_objectFromFile("OSMAPOSLReconstruction", "");
		handle = charDataHandle("reconstructedImage");
		RUN(cSTIR_setParameter
			(recon, "Reconstruction", "output_filename_prefix", handle));
		deleteDataHandle(handle);
		handle = intDataHandle(12);
		RUN(cSTIR_setParameter
			(recon, "IterativeReconstruction", "num_subsets", handle));
		deleteDataHandle(handle);
		handle = intDataHandle(num_subiterations);
		RUN(cSTIR_setParameter
			(recon, "IterativeReconstruction", "num_subiterations", handle));
		RUN(cSTIR_setParameter
			(recon, "IterativeReconstruction", "save_interval", handle));
		deleteDataHandle(handle);
		handle = intDataHandle(1);
		RUN(cSTIR_setParameter
			(recon, "IterativeReconstruction", "inter_iteration_filter_interval", 
			handle));
		deleteDataHandle(handle);
		RUN(cSTIR_setParameter
			(recon, "IterativeReconstruction", "objective_function", obj_fun));
		RUN(cSTIR_setParameter
			(recon, "IterativeReconstruction", "inter_iteration_filter_type", filter));
		handle = charDataHandle("multiplicative");
		RUN(cSTIR_setParameter(recon, "OSMAPOSL", "MAP_model", handle));
		deleteDataHandle(handle);
		RUN(cSTIR_setupReconstruction(recon, image));

		void* new_fd = cSTIR_acquisitionModelFwd(am, image, "");
		deleteDataHandle(new_fd);

		//for (int iter = 0; iter < num_subiterations; iter++) {
		//	std::cout << "iteration " << iter << '\n';
		//	cSTIR_updateReconstruction(recon, image);
		//}

		//double* rec_image_data = new double[image_size];
		//cSTIR_getImageData(image, (size_t)rec_image_data);
		//std::cout << "images diff: " << diff(image_size, image_data, rec_image_data)
		//	<< '\n';

		//delete[] rec_image_data;
		delete[] image_data;

		break;
	}

	deleteDataHandle(image);
	deleteDataHandle(matrix);
	deleteDataHandle(am);
	deleteDataHandle(ad);
	deleteDataHandle(fd);
	deleteDataHandle(at);
	deleteDataHandle(bt);
	deleteDataHandle(nd);
	deleteDataHandle(prior);
	deleteDataHandle(obj_fun);
	deleteDataHandle(filter);
	deleteDataHandle(recon);
}

