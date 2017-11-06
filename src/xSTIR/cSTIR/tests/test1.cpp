#include <string>

#include "object.h"
#include "stir_types.h"
#include "stir_x.h"
#include "SIRF/common/envar.h"

int test1()
{
	try {
		TextWriter w;
		TextWriterHandle h;
		h.set_information_channel(&w);

		std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

		std::string path(SIRF_path + "/examples/Python/PET");
		std::string filename;
		int dim[3];
		size_t sinos, views, tangs;

		//filename = SIRF_path + "/examples/Python/PET/my_image.hv";
		filename = SIRF_path + "/examples/Python/PET/reconstructed_image_2.hv";
		CREATE_OBJ(PETImageData, image_data, sptr_id, filename);

		image_data.get_dimensions(dim);
		size_t image_size = dim[0] * dim[1] * dim[2];
		std::cout << dim[0] << '\n';
		std::cout << dim[1] << '\n';
		std::cout << dim[2] << '\n';
		float im_norm = image_data.norm();
		std::cout << "image norm: " << im_norm << '\n';

		//filename = SIRF_path + "/data/examples/PET/Utahscat600k_ca_seg4.hs";
		filename = SIRF_path + "/data/examples/PET/my_forward_projection.hs";
		CREATE_OBJECT(PETAcquisitionData, PETAcquisitionDataInFile, acq_data, sptr_ad,
			filename.c_str());
		sinos = acq_data.get_num_sinograms();
		views = acq_data.get_num_views();
		tangs = acq_data.get_num_tangential_poss();
		float acq_norm = acq_data.norm();
		std::cout << "sinograms: " << sinos << '\n';
		std::cout << "views: " << views << '\n';
		std::cout << "tangential positions: " << tangs << '\n';
		std::cout << "acquisition data norm: " << acq_norm << '\n';

		shared_ptr<PETAcquisitionData> sptr_a = acq_data.new_acquisition_data();
		PETAcquisitionData& at = *sptr_a;
		at.fill(0.05f);
		shared_ptr<PETAcquisitionData> sptr_b = acq_data.new_acquisition_data();
		PETAcquisitionData& bt = *sptr_b;
		bt.fill(0.1f);
		shared_ptr<PETAcquisitionData> sptr_e = acq_data.new_acquisition_data();
		PETAcquisitionData& be = *sptr_e;
		be.fill(2.0f);

		CREATE_OBJ(RayTracingMatrix, matrix, sptr_matrix,);
		matrix.set_num_tangential_LORs(2);

		CREATE_OBJECT(PETAcquisitionModel, PETAcquisitionModelUsingMatrix, am, sptr_am,);
		am.set_matrix(sptr_matrix);
		//am.set_additive_term(sptr_a);
		//am.set_background_term(sptr_b);
		//am.set_bin_efficiency(sptr_e);
		am.set_up(sptr_ad, sptr_id);

		std::cout << "projecting...\n";
		shared_ptr<PETAcquisitionData> sptr_sd = am.forward(image_data);
		PETAcquisitionData& sim_data = *sptr_sd;
		sinos = sim_data.get_num_sinograms();
		views = sim_data.get_num_views();
		tangs = sim_data.get_num_tangential_poss();
		float sim_norm = sim_data.norm();
		std::cout << "sinograms: " << sinos << '\n';
		std::cout << "views: " << views << '\n';
		std::cout << "tangential positions: " << tangs << '\n';
		std::cout << "simulated data norm: " << sim_norm << '\n';

		//shared_ptr<aDataContainer<float> > sptr_z(acq_data.new_data_container());
		shared_ptr<PETAcquisitionData> sptr_diff(acq_data.new_acquisition_data());
		PETAcquisitionData& acq_diff = *sptr_diff;
		acq_diff.axpby(1.0/sim_norm, sim_data, -1.0/acq_norm, acq_data);
		std::cout << "relative difference: " << acq_diff.norm() << std::endl;

		shared_ptr<PETImageData> sptr_bd = am.backward(sim_data);
		PETImageData& back_data = *sptr_bd;
		back_data.get_dimensions(dim);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		float bd_norm = back_data.norm();
		std::cout << "back projection norm: " << bd_norm << '\n';
		shared_ptr<PETImageData> sptr_imd(image_data.new_image_data());
		PETImageData& img_diff = *sptr_imd;
		img_diff.axpby(1.0 / im_norm, image_data, -1.0 / bd_norm, back_data);
		std::cout << "relative difference: " << img_diff.norm() << std::endl;

		CREATE_OBJECT(Prior3DF, QuadPrior3DF, prior, sptr_prior,);
		prior.set_penalisation_factor(0.00f);

		CREATE_OBJECT(DataProcessor3DF, CylindricFilter3DF, filter, sptr_filter,);

		CREATE_OBJECT(ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF, obj_fun,
			sptr_fun,);
		obj_fun.set_acquisition_data(sptr_ad);
		obj_fun.set_acquisition_model(sptr_am);
		obj_fun.set_zero_seg0_end_planes(true);
		obj_fun.set_max_segment_num_to_process(3);
		obj_fun.set_prior_sptr(sptr_prior);

		int num_subiterations = 2;
		CREATE_OBJECT(Reconstruction<Image3DF>, OSMAPOSLReconstruction<Image3DF>, recon,
			sptr_recon,);
		recon.set_MAP_model("multiplicative");
		recon.set_num_subsets(12);
		recon.set_num_subiterations(num_subiterations);
		recon.set_save_interval(num_subiterations);
		recon.set_inter_iteration_filter_interval(1);
		recon.set_output_filename_prefix("reconstructedImage");
		recon.set_objective_function_sptr(sptr_fun);
		recon.set_inter_iteration_filter_ptr(sptr_filter);

		////Succeeded s = 
		////	//recon.set_up(sptr_image);
		////	xSTIR_setupReconstruction((void*)&sptr_recon, sptr_image);
		////if (s != Succeeded::yes) {
		////	std::cout << "xSTIR_setupReconstruction failed\n";
		////}
		xSTIR_IterativeReconstruction3DF& xrecon =
			(xSTIR_IterativeReconstruction3DF&)(recon);
		////(xSTIR_IterativeReconstruction3DF&)(*sptr_recon);
		Succeeded s = Succeeded::no;
		if (!xrecon.post_process()) {
			s = xrecon.setup(sptr_id->data_sptr());
			xrecon.subiteration() = xrecon.get_start_subiteration_num();
		}

		for (int iter = 0; iter < num_subiterations; iter++) {
			std::cout << "iteration " << iter << '\n';
			xrecon.update(sptr_id->data());
			//xrecon.update(image);
			//xSTIR_updateReconstruction((void*)&sptr_recon, image);
		}

		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		//double* rec_image_data = new double[image_size];
		//xSTIR_getImageDataAsDoubleArray(image, rec_image_data);
		//std::cout << "images diff: " << diff(image_size, image_data, rec_image_data) << '\n';
		//delete[] rec_image_data;
		//delete[] image_data;
	}
	catch (...) {
		std::cout << "exception thrown\n";
	}
	return 0;
}