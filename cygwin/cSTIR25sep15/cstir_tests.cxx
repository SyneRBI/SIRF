#include <iostream>
using namespace std;

#include "cstir.h"

void cstir_test1() {

	int not_found;
	void* handle;

	void* h_mx = cSTIR_newObject("RayTracingMatrix");
	not_found = cSTIR_setParameter
		(h_mx, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2));
	if (not_found)
		std::cout << "parameter not found" << std::endl;

	void* h_proj = cSTIR_newObject("ProjectorsUsingMatrix");
	not_found = cSTIR_setParameter
		(h_proj, "ProjectorsUsingMatrix", "matrix_type", h_mx);
	if (not_found)
		std::cout << "parameter not found" << std::endl;

	void* h_prior = cSTIR_newObject("QuadraticPrior");
	not_found = cSTIR_setParameter
		(h_prior, "GeneralisedPrior", "penalisation_factor", floatDataHandle(0.5));
	handle = cSTIR_parameter(h_prior, "GeneralisedPrior", "penalisation_factor");
	//std::cout << floatDataFromHandle(handle) << std::endl;
	deleteDataHandle(handle);

	void* h_filter = cSTIR_newObject("TruncateToCylindricalFOVImageProcessor");

	void* h_obj = cSTIR_newObject
		("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
	if (executionStatus(h_obj))
		std::cout << executionError(h_obj) << std::endl;
	not_found = cSTIR_setParameter
		(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
		"input_filename", charDataHandle("Utahscat600k_ca_seg4.hs"));
	if (not_found)
		std::cout << "parameter not found" << std::endl;
	not_found = cSTIR_setParameter
		(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
		"sensitivity_filename", charDataHandle("RPTsens_seg3_PM.hv"));
	not_found = cSTIR_setParameter
		(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
		"use_subset_sensitivities", charDataHandle("false"));
	not_found = cSTIR_setParameter
		(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
		"zero_seg0_end_planes", charDataHandle("true"));
	not_found = cSTIR_setParameter
		(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
		"max_segment_num_to_process", intDataHandle(3));
	not_found = cSTIR_setParameter
		(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
		"projector_pair_type", h_proj);
	not_found = cSTIR_setParameter
		(h_obj, "GeneralisedObjectiveFunction", "prior", h_prior);

	cSTIR_setupObjective(h_obj);

	void* h_recon = cSTIR_newReconstruction("OSMAPOSL", "");
	if (executionStatus(h_recon))
		std::cout << executionError(h_recon) << std::endl;

	not_found = cSTIR_setParameter
		(h_recon, "Reconstruction", "output_filename_prefix",
		charDataHandle("reconstructedImage"));

	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "num_subsets", intDataHandle(12));
	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "num_subiterations", intDataHandle(6));
	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "save_interval", intDataHandle(6));
	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "inter_iteration_filter_interval",
		intDataHandle(1));
	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "inter_iteration_filter_type",
		h_filter);
	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "objective_function", h_obj);
	not_found = cSTIR_setParameter
		(h_recon, "IterativeReconstruction", "initial_estimate",
		charDataHandle("my_uniform_image_circular.hv"));
	if (not_found)
		std::cout << "parameter not found" << std::endl;

	not_found = cSTIR_setParameter
		(h_recon, "OSMAPOSL", "MAP_model", charDataHandle("multiplicative"));

	void* h_image = cSTIR_imageFromFile("my_uniform_image_circular.hv");
	void* h_ximage = cSTIR_imageFromFile("test_image_PM_QP_6.hv");

	//handle = cSTIR_setupReconstruction("OSMAPOSL", h_recon, h_image);
	handle = cSTIR_setupReconstruction(h_recon, h_image);
	if (executionStatus(handle))
		std::cout << executionError(handle) << std::endl;
	deleteDataHandle(handle);
	handle = cSTIR_reconstruct(h_recon, h_image);
	if (executionStatus(handle))
		std::cout << executionError(handle) << std::endl;
	deleteDataHandle(handle);
	handle = cSTIR_imagesDifference(h_ximage, h_image, -1);
	if (executionStatus(handle))
		std::cout << executionError(handle) << std::endl;
	double diff = doubleDataFromHandle(handle);
	std::cout << "images difference: " << diff << std::endl;
	deleteDataHandle(handle);

	deleteDataHandle(h_mx);
	deleteDataHandle(h_proj);
	deleteDataHandle(h_prior);
	deleteDataHandle(h_filter);
	deleteDataHandle(h_obj);
	deleteDataHandle(h_recon);
	deleteDataHandle(h_image);
}
