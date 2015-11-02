//#include "stir/OSMAPOSL/OSMAPOSLReconstruction.h"
//#include "stir/DiscretisedDensity.h"
//#include "stir/Succeeded.h"
//#include "stir/CPUTimer.h"
//#include "stir/HighResWallClockTimer.h"
//#include "stir/recon_buildblock/distributable_main.h"

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

//extern "C"
//int doOSMAPOSLReconstruction(char* parFile);

//	USING_NAMESPACE_STIR

//	try {

//OSMAPOSLReconstruction<DiscretisedDensity<3, float> >
//	reconstruction_object(argc > 1 ? argv[1] : "OSMAPOSL_test_PM_QP.par0");

//if (reconstruction_object.reconstruct() == Succeeded::yes) {
//	_getch();
//	return EXIT_SUCCESS;
//}
//else {
//	_getch();
//	return EXIT_FAILURE;
//}
//}
//catch (...) {
//	cout << "exception occured" << endl;
//	_getch();
//}

//extern "C"
//int doOSMAPOSLReconstruction(char* parFile) {
//
//	USING_NAMESPACE_STIR
//
//	try {
//		OSMAPOSLReconstruction<DiscretisedDensity<3, float> >
//			reconstruction_object(parFile);
//		if (reconstruction_object.reconstruct() == Succeeded::yes) {
//			return EXIT_SUCCESS;
//		}
//		else {
//			return EXIT_FAILURE;
//		}
//	}
//	catch (...) {
//		cout << "exception occured" << endl;
//		return -1;
//	}
//}

//#include "dh.h"
//#include "ci_dh.h"
//#include "ci_ex.h"
//#include "ci_rp.h"
//#include "ci_stir.h"
//#include "ci_tw.h"
//#include "nlist.h"

//void* image;
//void* old_image;
//void* new_image;
//void* result;
//void* par;
//double diff;

////dh = runOSMAPOSLReconstruction(argc > 1 ? argv[1] : "OSMAPOSL_test_PM_QP.par");
////dh = STIR_OSMAPOSLReconstructionResult(argc > 1 ? argv[1] : "OSMAPOSL_test_PM_QP.par");
//image = STIR_OSMAPOSLReconstruction("OSMAPOSL_test_PM_QP.par");
//cout << executionStatus(image) << endl;
//cout << executionError(image) << endl;
//cout << executionErrorFile(image) << endl;
//cout << executionErrorLine(image) << endl;

////image = STIR_imageFromFile("my_test_image_PM_QP_6.hv");
//old_image = STIR_imageFromFile("test_image_PM_QP_6.hv");
//result = STIR_imageComparisonResult(old_image, image, -1, 0.00001);
//int same = intDataFromHandle(result);
//cout << same << endl;

//STIR_deleteImage(image);
//STIR_deleteImage(old_image);
//deleteTextWriter(tw);
//deleteTextWriterHandle(twh);

//par = STIR_parametersFromFile("test.par");
//par = STIR_parametersFromFile("OSMAPOSL_test_PM_QP.par");
//par = STIR_parametersFromFile("OSMAPOSL_test_PM_MRP.par");
//par = STIR_parametersFromFile("OSSPS_test_PM_QP.par");
par = STIR_parametersFromFile("OSMAPOSL_test.par");
if (executionStatus(par)) {
	cerr << executionError(par) << endl;
	cerr << "detected at line " << executionErrorLine(par);
	cerr << " of " << executionErrorFile(par) << endl;
}
else {
	STIR_printParameters(par);
	cout << STIR_parameter(par, "objective function type.input file") << endl;
	//STIR_printParameters(((DataHandle*)result)->data());

	//par = STIR_OSMAPOSLReconstructionParameters();
	//STIR_printParameters(par);

	image = STIR_reconstructedImage("OSMAPOSL", par);
	int dim[3];
	STIR_getImageDimensions(image, (size_t)dim);
	int nz = dim[0];
	int ny = dim[1];
	int nx = dim[2];
	cout << "nx: " << nx << endl;
	cout << "ny: " << ny << endl;
	cout << "nz: " << nz << endl;
	//double* data = STIR_imageData(image);
	double* data = new double[nx*ny*nz];
	STIR_getImageData(image, (size_t)data);
	//STIR_displayImageData(image);
	//image = STIR_OSMAPOSLReconstruction("OSMAPOSL_test.par");
	//image = STIR_OSMAPOSLReconstruction("OSMAPOSL_test_PM_MRP.par");
	//image = STIR_OSSPSReconstruction("OSSPS_test_PM_QP.par");
	//image = STIR_OSSPSReconstructedImage(par);
	//old_image = STIR_imageFromFile("test_image_OSSPS_PM_QP_8.hv");
	//old_image = STIR_imageFromFile("test_image_PM_MRP_6.hv");
	old_image = STIR_imageFromFile("test_image_5.hv");
	new_image = STIR_imageFromFile("my_test_image_5.hv");
	result = STIR_imagesDifference(old_image, new_image, -1);
	diff = doubleDataFromHandle(result);
	cout << "images difference: " << diff << endl;
	result = STIR_imagesDifference(image, new_image, -1);
	diff = doubleDataFromHandle(result);
	cout << "images difference: " << diff << endl;

	delete[] data;
	STIR_deleteImage(image);
	STIR_deleteImage(old_image);
	STIR_deleteImage(new_image);
}

STIR_deleteParameters(par);

//time_t seconds;
//seconds = time(NULL);
//char buff[64];
//sprintf_s(buff, 64, "tmp%d.par", seconds);
//cout << buff << endl;

//	PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> {

//#include "ci.h"
//#include "cstir.h"

//void* twh = newTextWriterHandle();
//void* tw = newTextWriter("try_stir.txt");
//setWriter(twh, tw);

xSTIR_Reconstruction(const char* parFile) : R(parFile) {}
void update(Image3DF &image) {
	update_estimate(image);
	end_of_iteration_processing(image);
	subiteration_num++;
}
int& subiteration() {
	return subiteration_num;
}
int subiteration() const {
	return subiteration_num;
}

//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>> sptr_obj1;
//sptr_obj1 = recon1.get_objective_function_sptr();
//obj = (xSTIR_theObjectiveFunction*)sptr_obj1.get();
//obj->set_input_file("Utahscat600k_ca_seg4.hs");
//obj->set_sensitivity_filename("RPTsens_seg3_PM.hv");
//obj->set_use_subset_sensitivities(false);
//obj->set_zero_seg0_end_planes(true);
//obj->set_max_segment_num_to_process(3);
//obj->set_projector_pair_sptr(sptr_proj);
//obj->set_prior_sptr(sptr_prior);

//xOSMAPOSLReconstruction recon1("OSMAPOSL_test_PM_QP.par");
//xOSMAPOSLReconstruction recon1("template.par");

//recon.set_objective_function_sptr(sptr_obj1);

//DataHandle* handle = (DataHandle*)h_obj;
//ptr = handle->data();
//boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>>* sptr_obj =
//	(boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>>*)ptr;
//xSTIR_theObjectiveFunction* obj = (xSTIR_theObjectiveFunction*)sptr_obj->get();
//obj->post_process();

//typedef xSTIR_ObjectiveFunction<PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>>
//typedef xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF
//xSTIR_theObjectiveFunction;

//handle = cSTIR_setParameter
//	(h_mx, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(-2));
//if (executionStatus(handle)) {
//	std::cout << executionError(handle) << std::endl;
//	deleteDataHandle(handle);
//	deleteDataHandle(h_mx);
//	return;
//}
//deleteDataHandle(handle);

//status = cSTIR_setParameter
//	(h_proj, "ProjectorsUsingMatrix", "matrix_type", h_mx);
////if (status)
////	std::cout << "parameter not found" << std::endl;

//status = cSTIR_setParameter
//	(h_prior, "GeneralisedPrior", "penalisation_factor", floatDataHandle(0.5));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_prior, "QuadraticPrior", "only_2D", intDataHandle(0));
////if (status)
////	std::cout << "parameter not found" << std::endl;
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "GeneralisedObjectiveFunction", "prior", h_prior);
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
//	"sensitivity_filename", charDataHandle("RPTsens_seg3_PM.hv"));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
//	"use_subset_sensitivities", charDataHandle("false"));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"input_filename", charDataHandle("Utahscat600k_ca_seg4.hs"));
////if (status)
////	std::cout << "parameter not found" << std::endl;
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"zero_seg0_end_planes", charDataHandle("true"));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"max_segment_num_to_process", intDataHandle(3));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"projector_pair_type", h_proj);
//deleteDataHandle(status);

//cSTIR_setupObjective(h_obj);

//status = cSTIR_setParameter
//	(h_recon, "Reconstruction", "output_filename_prefix",
//	charDataHandle("reconstructedImage"));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "num_subsets", intDataHandle(12));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "num_subiterations", intDataHandle(6));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "save_interval0", intDataHandle(6));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "inter_iteration_filter_interval",
//	intDataHandle(1));
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "inter_iteration_filter_type",
//	h_filter);
//deleteDataHandle(status);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "objective_function", h_obj);
//deleteDataHandle(status);
//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "initial_estimate",
//	charDataHandle("my_uniform_image_circular.hv"));
//if (status)
//	std::cout << "parameter not found" << std::endl;

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "inter_iteration_filter_interval",
//	intDataHandle(1));
//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "inter_iteration_filter_type",
//	h_filter);

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "initial_estimate",
//	charDataHandle("my_uniform_image_circular.hv"));
//if (status)
//	std::cout << "parameter not found" << std::endl;

//		("PoissonLogLikelihoodWithLinearModelForMeanAndProjData", h_obj);

//if (status)
//	std::cout << "parameter not found" << std::endl;

//status = cSTIR_setParameter
//	(h_recon, "OSMAPOSL", "MAP_model", charDataHandle("multiplicative"));
//deleteDataHandle(status);

//h_mx = cSTIR_newObject("RayTracingMatrix");
//status = cSTIR_setParameter
//	(h_mx, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2));

//h_proj = cSTIR_newObject("ProjectorsUsingMatrix");
//status = cSTIR_setParameter
//	(h_proj, "ProjectorsUsingMatrix", "matrix_type", h_mx);

//h_prior = cSTIR_newObject("QuadraticPrior");
//status = cSTIR_setParameter
//	(h_prior, "GeneralisedPrior", "penalisation_factor", floatDataHandle(0.5));
//handle = cSTIR_parameter(h_prior, "GeneralisedPrior", "penalisation_factor");
//std::cout << floatDataFromHandle(handle) << std::endl;
//deleteDataHandle(handle);

//h_filter = cSTIR_newObject("TruncateToCylindricalFOVImageProcessor");

//h_obj = cSTIR_newObject
//	("PoissonLogLikelihoodWithLinearModelForMeanAndProjData");
//if (executionStatus(h_obj))
//	std::cout << executionError(h_obj) << std::endl;
//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"input_filename", charDataHandle("Utahscat600k_ca_seg4.hs"));
//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
//	"sensitivity_filename", charDataHandle("RPTsens_seg3_PM.hv"));
//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMean",
//	"use_subset_sensitivities", charDataHandle("false"));
//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"zero_seg0_end_planes", charDataHandle("true"));
//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"max_segment_num_to_process", intDataHandle(3));
//status = cSTIR_setParameter
//	(h_obj, "PoissonLogLikelihoodWithLinearModelForMeanAndProjData",
//	"projector_pair_type", h_proj);
//status = cSTIR_setParameter
//	(h_obj, "GeneralisedObjectiveFunction", "prior", h_prior);

//cSTIR_setupObjective(h_obj);

//status = cSTIR_setParameter
//	(h_recon, "Reconstruction", "output_filename_prefix",
//	charDataHandle("reconstructedImage"));

//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "num_subsets", intDataHandle(4));
//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "num_subiterations", intDataHandle(8));
//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "save_interval", intDataHandle(8));
//status = cSTIR_setParameter
//	(h_recon, "IterativeReconstruction", "objective_function", h_obj);

//status = cSTIR_setParameter
//	(h_recon, "OSSPS", "relaxation_parameter", floatDataHandle(2.0));

//handle = cSTIR_setupReconstruction(h_recon, h_image);
//if (executionStatus(handle))
//	std::cout << executionError(handle) << std::endl;
//deleteDataHandle(handle);
//handle = cSTIR_reconstruct(h_recon, h_image);
//if (executionStatus(handle))
//	std::cout << executionError(handle) << std::endl;
//deleteDataHandle(handle);
//handle = cSTIR_imagesDifference(h_ximage, h_image, -1);
//if (executionStatus(handle))
//	std::cout << executionError(handle) << std::endl;
//double diff = doubleDataFromHandle(handle);
//std::cout << "images difference: " << diff << std::endl;
//deleteDataHandle(handle);

//deleteDataHandle(h_mx);
//deleteDataHandle(h_proj);
//deleteDataHandle(h_prior);
//deleteDataHandle(h_filter);
//deleteDataHandle(h_obj);
//deleteDataHandle(h_recon);
//deleteDataHandle(h_image);
//deleteDataHandle(h_ximage);

//deleteDataHandle(handle);
//handle = 

//if (!information_channel_)
//	std::cout << "set_information_channel: no information channel\n";

//if (!information_channel_)
//	std::cout << "set_warning_channel: no information channel\n";

//std::cout << "in print_information\n";
//if (!information_channel_)
//	std::cout << "no information channel\n";

//else
//	std::cout << "no information channel\n";

//std::cout << "initialized = " << initialized << std::endl;

//void* newTextWriterHandle() {
//	return (void*)new TextWriterHandle;
//}

//std::cout << channel << std::endl;
//std::cout << INFORMATION_CHANNEL << std::endl;
//std::cout << WARNING_CHANNEL << std::endl;

//void writeText(void* ptr, char* text) {
//	TextWriterHandle* h = (TextWriterHandle*)ptr;
//	h->write(text);
//}
//void deleteTextWriterHandle(void* ptr) {
//	delete (TextWriterHandle*)ptr;
//}

