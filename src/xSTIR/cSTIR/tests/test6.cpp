/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2022 Rutherford Appleton Laboratory STFC
Copyright 2018 - 2022 University College London
Copyright 2022 Harvard Medical School and Massachusetts General Hospital

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

\author Nikos Efthimiou
\author Kris Thielemans
*/
#include <cmath>
#include <fstream>
#include <string>

#include "stir/Verbosity.h"
#include "stir/IO/read_from_file.h"
#include "stir/find_STIR_config.h"
#include "sirf/common/iequals.h"
#include "sirf/STIR/stir_x.h"
#include "sirf/common/getenv.h"
#include "sirf/common/utilities.h"

#include "object.h"

using namespace sirf;

int test6()
{
        std::string STIR_data_path = stir::get_STIR_examples_dir();
        if (STIR_data_path.length() < 1) {
            std::cout << "cannot find data" << std::endl;
            return 1;
        }
        const std::string data_path = append_path(STIR_data_path, "samples");
        const std::string f_listmode = append_path(data_path, "mMR_listmode.l.hdr");
        //const std::string mu_map_filename = data_path + "mu_map.hv";
        //const std::string f_template = append_path(data_path, "mMR_template_span11_small.hs");

        std::cout << "running test6.cpp (PET listmode recon) with files from "
            << data_path << "...\n";

	try {
                stir::TextWriter w; // create writer with no output
		TextWriterHandle h;
		//h.set_information_channel(&w); // suppress STIR info output

		const std::string cache_path = data_path;
		//std::string sens_filename = cache_path + "sens_0.hv";
#if 0
		std::string tmpl_projdata_filename =
		    data_path + "tmpl_scanner.hs";

		CREATE_OBJECT(PETAcquisitionData, PETAcquisitionDataInFile,
			      acq_data, sptr_ad, tmpl_projdata_filename.c_str());

		PETAcquisitionDataInMemory::set_as_template();
#endif
#if 0
		// create compatible image
		CREATE_OBJ(STIRImageData, image_data, sptr_id, mu_map_filename.c_str());
		int dim[10];
		image_data.get_dimensions(dim);
		size_t image_size = dim[0] * dim[1] * dim[2];
		std::cout << "Image dimensions: "
			<< dim[0] << 'x' << dim[1] << 'x' << dim[2] << '\n';
		image_data.fill(1.0);
#endif
		CREATE_OBJECT(ObjectiveFunction3DF,
			      PoissonLLhLinModMeanListDataProjMatBin3DF,
			      obj_fun, sptr_fun,);
                std::shared_ptr<stir::ListModeData> listmode_data_sptr =
                  stir::read_from_file<stir::ListModeData>(f_listmode);
                obj_fun.set_input_data(listmode_data_sptr);
                std::shared_ptr<stir::DiscretisedDensity<3,float>> target_sptr(obj_fun.construct_target_ptr());
                CREATE_OBJ(STIRImageData, image_data, sptr_id, target_sptr);
                Coord3DI new_size = {-1,-1,-1};
                Coord3DF zooms = {1.f,4.f,4.f};
                sptr_id->zoom_image(zooms, Coord3DF(0.F,0.F,0.F), new_size, stir::ZoomOptions::preserve_projections);


		//This will activate use of caching and therefore multi-threading
		std::cout << "Setting cache path..." << std::endl;
		obj_fun.set_cache_path(cache_path);
		std::cout << "Cache path: " << obj_fun.get_cache_path() << std::endl;
		std::cout << "Setting max cache size ..." << std::endl;
		obj_fun.set_cache_max_size(1500000);
		std::cout << "Max cache: " << obj_fun.get_cache_max_size() << std::endl;
                obj_fun.set_recompute_cache(true);
#if 0
                // this functionality was for skip_lm_input_file, but this is disabled for now
		obj_fun.set_skip_lm_input_file(true);
		std::cout << "Setting scanner template..." << std::endl;
		obj_fun.set_acquisition_data(sptr_ad);
#endif
		std::cout << "Setting max segment num to process. ..." << std::endl;
		//obj_fun.set_max_segment_num_to_process(1); // small number for faster test
		//obj_fun.set_skip_balanced_subsets(true); // set for faster test

		//std::cout << "Sensitivity images: " << obj_fun.get_subsensitivity_filenames() << std::endl;
                //obj_fun.set_use_subset_sensitivities(false);
                // obj_fun.set_sensitivity_filename("1");

		int num_subsets = 1;

		xSTIR_OSMAPOSLReconstruction3DF recon;
		recon.set_num_subsets(num_subsets);
		recon.set_num_subiterations(num_subsets);
		std::cout << "Setting objective function ..." << std::endl;

		recon.set_objective_function_sptr(sptr_fun);
		recon.set_save_interval(1);
		recon.set_output_filename_prefix(data_path + "/my_output");
		std::cout <<  "Setting up the reconstruction, please wait ..." << std::endl;

		Succeeded s = recon.set_up(sptr_id->data_sptr());
		if (s == Succeeded::no)
			return 1;
		std::cout << "Reconstructor set up." << std::endl;

		recon.subiteration() = recon.get_start_subiteration_num();
		recon.reconstruct(sptr_id->data_sptr());

		return 0;
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
