#include <fstream>
#include <string>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"
//USING_NAMESPACE_STIR
//USING_NAMESPACE_ECAT

#include "object.h"
//#include "stir_types.h"
#include "stir_x.h"
//#include "SIRF/common/envar.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

bool file_exists(std::string filename)
{
	std::ifstream file;
	file.open(filename.c_str());
	if (file.good()) {
		file.close();
		return true;
	}
	return false;
}

int test1()
{
	try {
		TextWriter w;
		TextWriterHandle h;
		h.set_information_channel(&w);

		//std::string SIRF_path = EnvironmentVariable("SIRF_PATH");
		std::string SIRF_path = std::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

		std::string filename;
		int dim[3];
		size_t sinos, views, tangs;

		// if this file exists, the image stored there  will be used for comparison
		// with the reconstructed one
		filename = "reconstructedImage_2.hv";
		PETImageData* ptr_ei = 0;
		if (file_exists(filename)) {
			ptr_ei = new PETImageData(filename);
			ptr_ei->get_dimensions(dim);
			int nx = dim[2];
			int ny = dim[1];
			int nz = dim[0];
			std::cout << "image dimensions: "
				<< nx << 'x' << ny << 'x' << nz << '\n';
			PETImageData& image = *ptr_ei;
			float* ptr_data = new float[nx*ny*nz];
			image.get_data(ptr_data);
			int k = nx*ny*nz / 2;
			for (int i = k; i < k + 4; i++)
				std::cout << ptr_data[i] << '\n';
			auto iter = image.data().begin_all();
			for (int i = 0; i < k; i++, iter++);
			for (int i = k; i < k + 4; i++, iter++)
				std::cout << *iter << '\n';

		}
		return 0;
		// locate acquisition data
		//filename = SIRF_path + "/data/examples/PET/Utahscat600k_ca_seg4.hs";
		filename = SIRF_path + "/data/examples/PET/my_forward_projection.hs";
		CREATE_OBJECT(PETAcquisitionData, PETAcquisitionDataInFile, 
			acq_data, sptr_ad, filename.c_str());
		sinos = acq_data.get_num_sinograms();
		views = acq_data.get_num_views();
		tangs = acq_data.get_num_tangential_poss();
		float acq_norm = acq_data.norm();
		std::cout << "sinograms: " << sinos << '\n';
		std::cout << "views: " << views << '\n';
		std::cout << "tangential positions: " << tangs << '\n';
		std::cout << "acquisition data norm: " << acq_norm << '\n';

		// all acquisition-like data except acq_data will be stored in memory
		// (default storage is temporary files)
		PETAcquisitionDataInMemory::set_as_template();

		// create compatible image
		CREATE_OBJ(PETImageData, image_data, sptr_id, acq_data);
		image_data.get_dimensions(dim);
		size_t image_size = dim[0] * dim[1] * dim[2];
		std::cout << "image dimensions: "
			<< dim[0] << 'x' << dim[1] << 'x' << dim[2] << '\n';
		image_data.fill(1.0);
		float im_norm = image_data.norm();
		std::cout << "image norm: " << im_norm << '\n';

		// create additive term
		shared_ptr<PETAcquisitionData> sptr_a = acq_data.new_acquisition_data();
		PETAcquisitionData& at = *sptr_a;
		at.fill(0.05f);
		// create background term
		shared_ptr<PETAcquisitionData> sptr_b = acq_data.new_acquisition_data();
		PETAcquisitionData& bt = *sptr_b;
		bt.fill(0.1f);
		// create bin efficiencies term
		shared_ptr<PETAcquisitionData> sptr_e = acq_data.new_acquisition_data();
		PETAcquisitionData& be = *sptr_e;
		be.fill(2.0f);

		// create ray tracing matrix
		CREATE_OBJ(RayTracingMatrix, matrix, sptr_matrix,);
		matrix.set_num_tangential_LORs(2);

		// create acquisition model that uses ray tracing matrix
		CREATE_OBJECT(PETAcquisitionModel, PETAcquisitionModelUsingMatrix, 
			am, sptr_am,);
		am.set_matrix(sptr_matrix);
		am.set_additive_term(sptr_a);
		am.set_background_term(sptr_b);
		//am.set_bin_efficiency(sptr_e);
		am.set_up(sptr_ad, sptr_id);

		CREATE_OBJECT(ImageDataProcessor, xSTIR_SeparableGaussianImageFilter, procesor, sptr_processor,);
		processor.set_fwhms(stir:make_coords(3.F,4.F,3.F));
		am.set_image_data_processor(sptr_processor);

		// create quadratic prior
		CREATE_OBJECT(Prior3DF, QuadPrior3DF, prior, sptr_prior,);
		prior.set_penalisation_factor(0.00f);

		// create cylindric filter
		CREATE_OBJECT(DataProcessor3DF, CylindricFilter3DF, filter, sptr_filter,);

		// create Poisson log likelihood objective function
		CREATE_OBJECT(ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF, 
			obj_fun, sptr_fun,);
		obj_fun.set_acquisition_data(sptr_ad);
		obj_fun.set_acquisition_model(sptr_am);
		obj_fun.set_zero_seg0_end_planes(true);
		obj_fun.set_max_segment_num_to_process(4); // < 4 causes forward projection crash
		obj_fun.set_prior_sptr(sptr_prior);

		// create OSMAPOSL reconstructor
		int num_subiterations = 2;
		CREATE_OBJECT(Reconstruction3DF, OSMAPOSLReconstruction3DF, recon,
			sptr_recon, );
		recon.set_MAP_model("multiplicative");
		recon.set_num_subsets(12);
		recon.set_num_subiterations(num_subiterations);
		recon.set_save_interval(num_subiterations);
		recon.set_inter_iteration_filter_interval(1);
		recon.set_output_filename_prefix("reconstructedImage");
		recon.set_objective_function_sptr(sptr_fun);
		recon.set_inter_iteration_filter_ptr(sptr_filter);

		Succeeded s = recon.set_up(sptr_id);

		// reconstruct
		for (int iter = 0; iter < num_subiterations; iter++) {
			std::cout << "iteration " << recon.get_subiteration_num() << '\n';
			recon.update(sptr_id);
		}

		// forward-project the image to simulate the acquisition process
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

		// compare the simulated acquisition data with raw acquisition data
		shared_ptr<PETAcquisitionData> sptr_diff(acq_data.new_acquisition_data());
		PETAcquisitionData& acq_diff = *sptr_diff;
		acq_diff.axpby
			(float(1.0 / sim_norm), sim_data, -float(1.0 / acq_norm), acq_data);
		std::cout << "relative data difference: " << acq_diff.norm() << std::endl;

		// backproject the simulated data
		shared_ptr<PETImageData> sptr_bd = am.backward(sim_data);
		PETImageData& back_data = *sptr_bd;
		back_data.get_dimensions(dim);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		float bd_norm = back_data.norm();
		std::cout << "back projection norm: " << bd_norm << '\n';

		// compare backprojected data with the reconstructed image
		shared_ptr<PETImageData> sptr_imd(image_data.new_image_data());
		PETImageData& img_diff = *sptr_imd;
		im_norm = image_data.norm();
		img_diff.axpby
			(float(1.0 / im_norm), image_data, -float(1.0 / bd_norm), back_data);
		std::cout << "relative images difference: " << img_diff.norm() << std::endl;

		if (ptr_ei) {
			// compare the reconstructed and expected images
			img_diff.axpby(1.0f, image_data, -1.0f, *ptr_ei);
			std::cout << "images difference: " << img_diff.norm() << std::endl;
		}
	}
	catch (...) {
		std::cout << "exception thrown\n";
	}
	return 0;
}
