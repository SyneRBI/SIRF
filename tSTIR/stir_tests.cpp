#include <iostream>
using namespace std;

#include "stir_x.h"

#define CAST_PTR(T, X, Y) T* X = (T*)Y

typedef CartesianCoordinate3D<float> Coord3DF;

//typedef xSTIR_ObjectiveFunction<PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF>>
typedef xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF
xSTIR_theObjectiveFunction;
//typedef xSTIR_IterativeReconstruction<OSMAPOSLReconstruction<Image3DF>>
//xOSMAPOSLReconstruction;

void stir_test0() {
	OSMAPOSLReconstruction<Image3DF> recon("OSMAPOSL_test_PM_QP.par");
	sptrImage3DF* sptr_image = new sptrImage3DF
		(read_from_file<Image3DF>("my_uniform_image_circular.hv"));
	xSTIR_setupReconstruction //<OSMAPOSLReconstruction<Image3DF>>
		(&recon, *sptr_image);
	//recon.setup(*sptr_image);
}

void stir_test1() {

	boost::shared_ptr<ProjMatrixByBin> sptr_mx(new ProjMatrixByBinUsingRayTracing);
	//sptr_mx->set_num_tangential_LORs(2);
	((ProjMatrixByBinUsingRayTracing*)sptr_mx.get())->set_num_tangential_LORs(2);

	boost::shared_ptr<ProjectorByBinPairUsingProjMatrixByBin> sptr_proj
		(new ProjectorByBinPairUsingProjMatrixByBin);
	sptr_proj->set_proj_matrix_sptr(sptr_mx);

	boost::shared_ptr<GeneralisedPrior<Image3DF>> sptr_prior(new QuadraticPrior<float>);
	((QuadraticPrior<float>*)sptr_prior.get())->set_penalisation_factor(0.5);

	boost::shared_ptr<DataProcessor<Image3DF>> sptr_filter
		(new TruncateToCylindricalFOVImageProcessor<float>);

	boost::shared_ptr<GeneralisedObjectiveFunction<Image3DF>> sptr_obj
		(new xSTIR_theObjectiveFunction);
	xSTIR_theObjectiveFunction* obj = (xSTIR_theObjectiveFunction*)sptr_obj.get();
	obj->set_input_file("Utahscat600k_ca_seg4.hs");
	obj->set_sensitivity_filename("RPTsens_seg3_PM.hv");
	//obj->set_recompute_sensitivity(true);
	obj->set_use_subset_sensitivities(false);
	obj->set_zero_seg0_end_planes(true);
	obj->set_max_segment_num_to_process(3);
	obj->set_projector_pair_sptr(sptr_proj);
	obj->set_prior_sptr(sptr_prior);
	((xSTIR_GeneralisedObjectiveFunction3DF*)obj)->post_process();

	OSMAPOSLReconstruction<Image3DF> recon;
	recon.set_MAP_model("multiplicative");
	recon.set_num_subsets(12);
	recon.set_num_subiterations(6);
	recon.set_save_interval(6);
	recon.set_inter_iteration_filter_interval(1);
	recon.set_inter_iteration_filter_ptr(sptr_filter);
	recon.set_output_filename_prefix("reconstructedImage");
	//recon.set_initial_estimate_file("my_uniform_image_circular.hv");
	xSTIR_set_initial_estimate_file //<IterativeReconstruction<Image3DF>>
		(&recon, "my_uniform_image_circular.hv");
	recon.set_objective_function_sptr(sptr_obj);

	//std::cout << "=======================================================" << std::endl;
	//std::cout << recon.parameter_info() << std::endl;;

	sptrImage3DF* sptr_image = new sptrImage3DF
		(read_from_file<Image3DF>("my_uniform_image_circular.hv"));
	sptrImage3DF* sptr_ximage = new sptrImage3DF
		(read_from_file<Image3DF>("test_image_PM_QP_6.hv"));

	xSTIR_setupReconstruction //<OSMAPOSLReconstruction<Image3DF>>
		(&recon, *sptr_image);
	//recon.setup(*sptr_image);
	recon.reconstruct(*sptr_image);

}

void stir_test2() {
	boost::shared_ptr<Shape3D> sptr(new EllipsoidalCylinder);
	Coord3DF my_origin(1, 2, 3);
	sptr->set_origin(my_origin);
	Coord3DF origin = sptr->get_origin();
	cout << "origin:" << endl;
	cout << origin.x() << endl;
	cout << origin.y() << endl;
	cout << origin.z() << endl;
	CAST_PTR(EllipsoidalCylinder, ptr, sptr.get());
	ptr->set_length(10);
	ptr->set_radius_x(20);
	ptr->set_radius_y(30);
	cout << "length:" << endl;
	cout << ptr->get_length() << endl;
	cout << "radii:" << endl;
	cout << ptr->get_radius_x() << endl;
	cout << ptr->get_radius_y() << endl;
}