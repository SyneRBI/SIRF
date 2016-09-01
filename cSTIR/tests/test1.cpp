#include <string>

#include "data_handle.h"
#include "stir.h"
#include "stir_x.h"
#include "tests.h"

void test1()
{
	try {
		OBJECT(ProjMatrixByBin, RayTracingMatrix, matrix, sptr_matrix);
		matrix.set_num_tangential_LORs(2);

		OBJECT(ProjectorByBinPair, ProjectorPairUsingMatrix, am, sptr_am);
		am.set_proj_matrix_sptr(sptr_matrix);

		std::string path("../../examples/");
		std::string filename(path + "my_forward_projection.hs");
		boost::shared_ptr<ProjData> sptr_ad = ProjData::read_from_file(filename);
		size_t size = sptr_ad->size_all();
		size_t sinos = sptr_ad->get_num_sinograms();
		size_t views = sptr_ad->get_num_views();
		size_t tang_pos = sptr_ad->get_num_tangential_poss();
		std::cout << "segments: " << sptr_ad->get_num_segments() << '\n';
		std::cout << "sinograms: " << sinos << '\n';
		std::cout << "views: " << views << '\n';
		std::cout << "tangential positions: " << tang_pos << '\n';
		std::cout << "size: " << size << ' ' << sinos*views*tang_pos << '\n';
		double* acq_data = new double[size];
		sptr_ad->copy_to(acq_data);
		delete[] acq_data;

		OBJECT(Prior3DF, QuadPrior3DF, prior, sptr_prior);
		prior.set_penalisation_factor(0.001f);

		OBJECT(DataProcessor3DF, CylindricFilter3DF, filter, sptr_filter);

		OBJECT(ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF, obj_fun,
			sptr_fun);
		obj_fun.set_zero_seg0_end_planes(true);
		obj_fun.set_max_segment_num_to_process(3);
		obj_fun.set_projector_pair_sptr(sptr_am);
		obj_fun.set_proj_data_sptr(sptr_ad);
		obj_fun.set_prior_sptr(sptr_prior);

		int num_subiterations = 6;
		OBJECT(Reconstruction<Image3DF>, OSMAPOSLReconstruction<Image3DF>, recon,
			sptr_recon);
		recon.set_MAP_model("multiplicative");
		recon.set_num_subsets(12);
		recon.set_num_subiterations(num_subiterations);
		recon.set_save_interval(num_subiterations);
		recon.set_inter_iteration_filter_interval(1);
		recon.set_output_filename_prefix("reconstructedImage");
		recon.set_objective_function_sptr(sptr_fun);
		recon.set_inter_iteration_filter_ptr(sptr_filter);

		filename = path + "expected_image.hv";
		sptrImage3DF sptr_image(read_from_file<Image3DF>(filename));
		Image3DF& image = *sptr_image;

		Succeeded s = xSTIR_setupReconstruction((void*)&sptr_recon, sptr_image);
		if (s != Succeeded::yes) {
			std::cout << "xSTIR_setupReconstruction failed\n";
		}

		for (int iter = 0; iter < num_subiterations; iter++) {
			std::cout << "iteration " << iter << '\n';
			xSTIR_updateReconstruction((void*)&sptr_recon, image);
		}
	}
	catch (...) {
		std::cout << "exception thrown\n";
	}
}