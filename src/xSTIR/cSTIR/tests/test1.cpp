/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 Rutherford Appleton Laboratory STFC
Copyright 2018 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

/*!
\file
\ingroup PET

\author Evgueni Ovtchinnikov
\author Richard Brown
\author SyneRBI
*/
#include <cmath>
#include <fstream>
#include <string>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"
#include "stir/ExamData.h"

#include "sirf/common/iequals.h"
#include "sirf/STIR/stir_x.h"

#include "getenv.h"
#include "object.h"

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

extern "C"
void openChannel(int channel, void* ptr_w);

int test1()
{
	std::cout << "running test1.cpp...\n";

	try {
		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

		TextWriter w; // create writer with no output
		TextWriterHandle h;
		h.set_information_channel(&w); // suppress STIR info output

		bool ok;
		bool fail = false;

		std::string filename;
		int dim[10];
		size_t sinos, views, tangs;
		// locate acquisition data
		//filename = SIRF_path + "/data/examples/PET/Utahscat600k_ca_seg4.hs";
		filename = SIRF_path + "/data/examples/PET/my_forward_projection.hs";
		fix_path_separator(filename);
		CREATE_OBJECT(PETAcquisitionData, PETAcquisitionDataInFile,
			acq_data, sptr_ad, filename.c_str());
		sinos = acq_data.get_num_sinograms();
		views = acq_data.get_num_views();
		//views = sptr_ad->get_num_views();
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
		CREATE_OBJ(STIRImageData, image_data, sptr_id, acq_data);
		image_data.get_dimensions(dim);
		size_t image_size = dim[0] * dim[1] * dim[2];
		std::cout << "image dimensions: "
			<< dim[0] << 'x' << dim[1] << 'x' << dim[2] << '\n';
		image_data.fill(1.0);
		float im_norm = image_data.norm();
		std::cout << "image norm: " << im_norm << '\n';
		const VoxelisedGeometricalInfo3D &geom_info = *image_data.get_geom_info_sptr();
		const VoxelisedGeometricalInfo3D &geom_info_copy = *image_data.get_geom_info_sptr();
		std::cout << geom_info.get_info().c_str();
		ok = (geom_info == geom_info_copy);
		if (ok)
			std::cout << "geom_info == ok\n";
		else
			std::cout << "geom_info == failed \n";
		fail = fail || !ok;

		// show and change modality demo
		std::string mod = image_data.modality();
		std::cout << '\n' << "modality: " << mod << '\n';
		image_data.set_modality("NM");
		std::cout << "new modality set: " << image_data.modality() << '\n';
		ok = sirf::iequals(image_data.modality(), "NM");
		if (ok)
			std::cout << "set_modality ok\n";
		else
			std::cout << "set_modality failed \n";
		fail = fail || !ok;
		// restore
		image_data.set_modality(mod);

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

		// create acquisition model that uses ray tracing matrix
		// long way:
		// create ray tracing matrix
		//CREATE_OBJ(RayTracingMatrix, matrix, sptr_matrix, );
		//matrix.set_num_tangential_LORs(2);
		//CREATE_OBJECT(PETAcquisitionModel, PETAcquisitionModelUsingMatrix,
		//	am, sptr_am, );
		//am.set_matrix(sptr_matrix);
		// short way:
		CREATE_OBJECT(PETAcquisitionModel, PETAcquisitionModelUsingRayTracingMatrix,
			am, sptr_am,);
		am.set_num_tangential_LORs(10);
		am.set_additive_term(sptr_a);
		am.set_background_term(sptr_b);
		am.set_up(sptr_ad, sptr_id);

		int num_LORs = am.get_num_tangential_LORs();
		std::cout << "tangential LORs: " << num_LORs << '\n';

		CREATE_OBJECT(ImageDataProcessor, xSTIR_SeparableGaussianImageFilter, processor, sptr_processor,);
//		processor.set_fwhms(stir::make_coords(3.F, 4.F, 3.F));
		stir::Coordinate3D< float > fwhms(3.F, 4.F, 3.F);
		processor.set_fwhms(fwhms);
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
		int num_subiterations = 4;
		xSTIR_OSMAPOSLReconstruction3DF recon;
		recon.set_MAP_model("multiplicative");
		recon.set_num_subsets(12);
		recon.set_num_subiterations(num_subiterations);
		recon.set_objective_function_sptr(sptr_fun);
		recon.set_input_data(sptr_ad->data());
		/* the stuff below is for recon.reconstruct() */
		//recon.set_save_interval(num_subiterations);
		//recon.set_inter_iteration_filter_interval(1);
		//recon.set_output_filename_prefix("reconstructedImage");
		//recon.set_inter_iteration_filter_ptr(sptr_filter);

		std::cout << "setting up the reconstructor, please wait...";
		Succeeded s = recon.set_up(sptr_id->data_sptr());
		std::cout << "ok\n";
		//recon.reconstruct(sptr_id->data_sptr());

		// reconstruct
		recon.subiteration() = recon.get_start_subiteration_num();
		for (int iter = 0; iter < num_subiterations; iter++) {
			std::cout << "iteration " << iter << '\n';
			recon.update_estimate(sptr_id->data());
			std::cout << "image norm: " << image_data.norm() << '\n';
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
		float alpha = 1.0 / sim_norm;
		float beta = -alpha;
		acq_diff.axpby
			(&alpha, sim_data, &beta, acq_data);
		std::cout << "relative data difference: " << acq_diff.norm() << std::endl;

		// backproject the simulated data
		shared_ptr<STIRImageData> sptr_bd = am.backward(sim_data);
		STIRImageData& back_data = *sptr_bd;
		back_data.get_dimensions(dim);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		float bd_norm = back_data.norm();
		std::cout << "back projection norm: " << bd_norm << '\n';

		// compare backprojected data with the reconstructed image
		shared_ptr<STIRImageData> sptr_imd(image_data.new_image_data());
		STIRImageData& img_diff = *sptr_imd;
		im_norm = image_data.norm();
		bd_norm = back_data.norm();
		alpha = 1.0 / im_norm;
		beta = -1.0 / bd_norm;
		img_diff.axpby
			(&alpha, image_data, &beta, back_data);
		std::cout << "relative images difference: " << img_diff.norm() << std::endl;

		// compute the norm of the linear part of the acquisition model
		std::cout << "\ncomputing the norm of the linear part of the acquisition model...\n";
		float am_norm = am.norm();

		std::cout << "\nchecking the acquisition model norm:\n";
		std::cout << "acquisition model norm: |A| = " << am_norm << '\n';
		std::cout << "image data x norm: |x| = " << im_norm << '\n';
		std::cout << "simulated acquisition data norm: |A(x)| = " << sim_norm << '\n';
		std::cout << "checking that |A(x)| <= |A||x|: ";
		float bound = am_norm*im_norm;
		ok = (sim_norm <= bound);
		if (ok)
			std::cout << sim_norm << " <= " << bound << " ok!\n";
		else
			std::cout << sim_norm << " > " << bound << " failure!\n";
		fail = fail || !ok;

		// restore the default storage scheme
		PETAcquisitionDataInFile::set_as_template();

		h.set_information_channel(0);

		std::cout << "done with test1.cpp...\n";

		return fail;
	}
	catch (const std::exception &error) {
		std::cerr << "\nException thrown:\n\t" << error.what() << "\n\n";
		return 1;
	}
	catch (...) {
		std::cerr << "\nException thrown\n";
		return 1;
	}
	return 0;
}
