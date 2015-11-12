#include <iostream>
using namespace std;

#include "dh.h"
#include "cstir.h"
#include "stir_x.h"

typedef CartesianCoordinate3D<float> Coord3DF;
typedef VoxelsOnCartesianGrid<float> Voxels3DF;

int setParameter(void* hs, const char* set, const char* par, void* hv)
{
	int status = 0;
	DataHandle* handle = (DataHandle*)cSTIR_setParameter(hs, set, par, hv);
	if (executionStatus(handle)) {
		std::cout << executionError(handle) << std::endl;
		status = 1;
	}
	deleteDataHandle(handle);
	return status;
}

void cstir_test1() {

	void* handle;

	void* h_mx = 0;
	void* h_proj = 0;
	void* h_prior = 0;
	void* h_filter = 0;
	void* h_obj = 0;
	void* h_recon = 0;
	void* h_image = 0;
	void* h_ximage = 0;

	void* cinf = newTextPrinter("stdout");
	void* cwrn = newTextWriter("wrn.txt");
	void* cerr = newTextWriter("err.txt");

	for (;;)  {
		h_mx = cSTIR_newObject("RayTracingMatrix");
		if (setParameter
			(h_mx, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2)))
			break;

		h_proj = cSTIR_newObject("ProjectorsUsingMatrix");
		if (setParameter(h_proj, "ProjectorsUsingMatrix", "matrix_type", h_mx))
			break;

		h_prior = cSTIR_newObject("QuadraticPrior");
		if (setParameter(h_prior, "GeneralisedPrior", "penalisation_factor", 
			floatDataHandle(0.5)))
			break;

		if (setParameter(h_prior, "QuadraticPrior", "only_2D", intDataHandle(0)))
			break;

		handle = cSTIR_parameter(h_prior, "GeneralisedPrior", "penalisation_factor");
		std::cout << floatDataFromHandle(handle) << std::endl;
		deleteDataHandle(handle);

		handle = cSTIR_setupObject("GeneralisedPrior", h_prior);
		if (executionStatus(handle)) {
			std::cout << executionError(handle) << std::endl;
			deleteDataHandle(handle);
			break;
		}
		deleteDataHandle(handle);

		h_filter = cSTIR_newObject("TruncateToCylindricalFOVImageProcessor");

		h_obj = cSTIR_newObject
			("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
		if (executionStatus(h_obj))
			std::cout << executionError(h_obj) << std::endl;

		if (setParameter(h_obj, "GeneralisedObjectiveFunction", "prior", h_prior))
			break;

		if (setParameter(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
			"sensitivity_filename", charDataHandle("RPTsens_seg3_PM.hv")))
			break;
		if (setParameter(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
			"use_subset_sensitivities", charDataHandle("false")))
			break;

		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"input_filename", charDataHandle("Utahscat600k_ca_seg4.hs")))
			break;
		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"zero_seg0_end_planes", charDataHandle("true")))
			break;
		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"max_segment_num_to_process", intDataHandle(3)))
			break;
		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"projector_pair_type", h_proj))
			break;

		handle = cSTIR_setupObject("GeneralisedObjectiveFunction", h_obj);
		if (executionStatus(handle)) {
			std::cout << executionError(handle) << std::endl;
			deleteDataHandle(handle);
			break;
		}
		deleteDataHandle(handle);

		//h_recon = cSTIR_newReconstruction("OSMAPOSL", "");
		h_recon = cSTIR_objectFromFile("OSMAPOSLReconstruction", "");
		if (executionStatus(h_recon))
			std::cout << executionError(h_recon) << std::endl;

		if (setParameter(h_recon, "Reconstruction", "output_filename_prefix",
			charDataHandle("reconstructedImage")))
			break;

		if (setParameter
			(h_recon, "IterativeReconstruction", "num_subsets", intDataHandle(12)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "start_subset_num", intDataHandle(0)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "num_subiterations", intDataHandle(6)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "save_interval", intDataHandle(6)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "inter_iteration_filter_interval",
			intDataHandle(1)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "inter_iteration_filter_type",
			h_filter))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "objective_function", h_obj))
			break;

		if (setParameter
			(h_recon, "OSMAPOSL", "MAP_model", charDataHandle("multiplicative")))
			break;

		h_image = cSTIR_objectFromFile("Image", "my_uniform_image_circular.hv");
		h_ximage = cSTIR_objectFromFile("Image", "test_image_PM_QP_6.hv");

		for (;;) {
			handle = cSTIR_setupReconstruction(h_recon, h_image);
			if (executionStatus(handle)){
				std::cout << executionError(handle) << std::endl;
				deleteDataHandle(handle);
				break;
			}
			deleteDataHandle(handle);
			handle = cSTIR_runReconstruction(h_recon, h_image);
			if (executionStatus(handle)) {
				std::cout << executionError(handle) << std::endl;
				deleteDataHandle(handle);
				break;
			}
			deleteDataHandle(handle);
			handle = cSTIR_imagesDifference(h_ximage, h_image, -1);
			if (executionStatus(handle)) {
				std::cout << executionError(handle) << std::endl;
				deleteDataHandle(handle);
				break;
			}
			double diff = doubleDataFromHandle(handle);
			std::cout << "images difference: " << diff << std::endl;
			deleteDataHandle(handle);

			char buff[64];
			setWriter(cinf, 0);
			setWriter(cwrn, 1);
			setWriter(cerr, 2);
			sprintf(buff, "images difference: %e\n", diff);
			writeText(buff);
			writeText(buff, WARNING_CHANNEL);
			writeText(buff, ERROR_CHANNEL);
			break;
		}
		break;
	}

	//cSTIR_deleteObject(h_mx, "RayTracingMatrix");
	//cSTIR_deleteObject(h_proj, "ProjectorsUsingMatrix");
	//cSTIR_deleteObject(h_prior, "GeneralisedPrior");
	//cSTIR_deleteObject(h_prior, "QuadraticPrior");
	cSTIR_deleteObject(h_filter); // , "DataProcessor");
	//	(h_filter, "TruncateToCylindricalFOVImageProcessor");
	//cSTIR_deleteObject
	//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
	cSTIR_deleteObject(h_mx); // , "ProjMatrix");
	cSTIR_deleteObject(h_proj); // , "Projectors");
	cSTIR_deleteObject(h_prior); // , "Prior");
	cSTIR_deleteObject(h_obj); // , "ObjectiveFunction");
	cSTIR_deleteObject(h_recon); // , "Reconstruction");
	cSTIR_deleteObject(h_image); // , "Image");
	cSTIR_deleteObject(h_ximage); // , "Image");

	deleteTextPrinter(cinf);
	deleteTextWriter(cwrn);
	deleteTextWriter(cerr);
}

void cstir_test2() {

	void* handle;

	void* h_mx = 0;
	void* h_proj = 0;
	void* h_prior = 0;
	void* h_obj = 0;
	void* h_recon = 0;
	void* h_image = 0;
	void* h_ximage = 0;

	for (;;) {

		h_mx = cSTIR_newObject("RayTracingMatrix");
		if (setParameter
			(h_mx, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2)))
			break;

		h_proj = cSTIR_newObject("ProjectorsUsingMatrix");
		if (setParameter(h_proj, "ProjectorsUsingMatrix", "matrix_type", h_mx))
			break;

		h_prior = cSTIR_newObject("QuadraticPrior");
		if (setParameter(h_prior, "GeneralisedPrior", "penalisation_factor",
			floatDataHandle(0.5)))
			break;

		if (setParameter(h_prior, "QuadraticPrior", "only_2D", intDataHandle(0)))
			break;

		handle = cSTIR_parameter(h_prior, "GeneralisedPrior", "penalisation_factor");
		std::cout << floatDataFromHandle(handle) << std::endl;
		deleteDataHandle(handle);

		handle = cSTIR_setupObject("prior", h_prior);
		if (executionStatus(handle))
			std::cout << executionError(handle) << std::endl;
		deleteDataHandle(handle);

		h_obj = cSTIR_newObject
			("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
		if (executionStatus(h_obj))
			std::cout << executionError(h_obj) << std::endl;

		if (setParameter(h_obj, "GeneralisedObjectiveFunction", "prior", h_prior))
			break;

		if (setParameter(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
			"sensitivity_filename", charDataHandle("RPTsens_seg3_PM.hv")))
			break;
		if (setParameter(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
			"use_subset_sensitivities", charDataHandle("false")))
			break;

		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"input_filename", charDataHandle("Utahscat600k_ca_seg4.hs")))
			break;
		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"zero_seg0_end_planes", charDataHandle("true")))
			break;
		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"max_segment_num_to_process", intDataHandle(3)))
			break;
		if (setParameter
			(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
			"projector_pair_type", h_proj))
			break;

		handle = cSTIR_setupObject("objective_function", h_obj);
		if (executionStatus(handle))
			std::cout << executionError(handle) << std::endl;
		deleteDataHandle(handle);

		//h_recon = cSTIR_newReconstruction("OSSPS", "");
		h_recon = cSTIR_objectFromFile("OSSPSReconstruction", "");
		if (executionStatus(h_recon))
			std::cout << executionError(h_recon) << std::endl;

		if (setParameter(h_recon, "Reconstruction", "output_filename_prefix",
			charDataHandle("reconstructedImage")))
			break;

		if (setParameter
			(h_recon, "IterativeReconstruction", "num_subsets", intDataHandle(4)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "num_subiterations", intDataHandle(8)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "save_interval", intDataHandle(8)))
			break;
		if (setParameter
			(h_recon, "IterativeReconstruction", "objective_function", h_obj))
			break;

		if (setParameter
			(h_recon, "OSSPS", "relaxation_parameter", floatDataHandle(2.0)))
			break;

		h_image = cSTIR_objectFromFile("Image", "test_image_PM_QP_6.hv");
		h_ximage = cSTIR_objectFromFile("Image", "test_image_OSSPS_PM_QP_8.hv");

		for (;;) {
			handle = cSTIR_setupReconstruction(h_recon, h_image);
			if (executionStatus(handle)){
				std::cout << executionError(handle) << std::endl;
				deleteDataHandle(handle);
				break;
			}
			deleteDataHandle(handle);
			handle = cSTIR_runReconstruction(h_recon, h_image);
			if (executionStatus(handle)) {
				std::cout << executionError(handle) << std::endl;
				deleteDataHandle(handle);
				break;
			}
			deleteDataHandle(handle);
			handle = cSTIR_imagesDifference(h_ximage, h_image, -1);
			if (executionStatus(handle)) {
				std::cout << executionError(handle) << std::endl;
				deleteDataHandle(handle);
				break;
			}
			double diff = doubleDataFromHandle(handle);
			std::cout << "images difference: " << diff << std::endl;
			deleteDataHandle(handle);
			break;
		}
		break;
	}

	cSTIR_deleteObject(h_mx); // , "ProjMatrix");
	cSTIR_deleteObject(h_proj); // , "Projectors");
	cSTIR_deleteObject(h_prior); // , "Prior");
	cSTIR_deleteObject(h_obj); // , "ObjectiveFunction");
	cSTIR_deleteObject(h_recon); // , "Reconstruction");
	cSTIR_deleteObject(h_image); // , "Image");
	cSTIR_deleteObject(h_ximage); // , "Image");
}

void cstir_test3() {

	void* h_v = 0;
	void* h_i = 0;
	void* h_s = 0;
	void* handle;
	int dim[3];

	for (;;) {
		h_v = cSTIR_voxels3DF(10, 20, 30, 1., 2., 3., 100., 200., 300.);
		h_i = cSTIR_imageFromVoxels(h_v);
		h_s = cSTIR_newObject("EllipsoidalCylinder");
		handle = cSTIR_addShape(h_i, h_v, h_s, 1.0);
		deleteDataHandle(handle);
		cSTIR_getImageDimensions(h_i, (size_t)dim);
		cout << dim[0] << endl;
		cout << dim[1] << endl;
		cout << dim[2] << endl;
		handle = cSTIR_setParameter(h_s, "Shape", "x", floatDataHandle(25));
		deleteDataHandle(handle);
		handle = cSTIR_parameter(h_s, "Shape", "x");
		cout << floatDataFromHandle(handle) << endl;
		deleteDataHandle(handle);
		break;
	}

	cSTIR_deleteObject(h_v); // , "Voxels");
	cSTIR_deleteObject(h_i); // , "Image");
	cSTIR_deleteObject(h_s); // , "Shape3D");
	//cSTIR_deleteObject(h_s, "EllipsoidalCylinder");
}