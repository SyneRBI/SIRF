/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2020 Rutherford Appleton Laboratory STFC
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
#include <iostream>
#include <cstdlib>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"
#include "stir/Verbosity.h"

#include "sirf/STIR/stir_x.h"
#include "sirf/common/getenv.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

int test4()
{
	std::cout << "running test4.cpp...\n";
	try {

		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

		std::string path = SIRF_path + "/data/examples/PET/mMR/";

		std::string f_listmode = path + "list.l.hdr";
		std::string f_template = path + "mMR_template_span11_small.hs";
		PETAcquisitionDataInFile acq_data_template(f_template.c_str());

		// Listmode to sinograms
		ListmodeToSinograms converter;
		converter.set_input(f_listmode);
		converter.set_output("proj_data");
		converter.set_template(acq_data_template);
		//// old way (now just an alternative option):
		//converter.set_template(f_template);
		converter.set_time_interval(0, 10);
		converter.set_up();
		converter.estimate_randoms();
		converter.save_randoms();

        // Check count rates - for the particular dataset,
        // we know that 73036 is exceeded at 22s. You can
        // see this with STIR's list_lm_countrates
        const unsigned long num_prompt_threshold = 73036.f;
        const float known_time = 22.f;
        const float time_at_which_num_prompts_exceeds_threshold =
                converter.get_time_at_which_num_prompts_exceeds_threshold(num_prompt_threshold);
        if (std::abs(time_at_which_num_prompts_exceeds_threshold-known_time) > 1e-4f)
            throw std::runtime_error("ListmodeToSinograms::get_time_at_which_num_prompts_exceeds_threshold failed");

        // Construct STIRImageData from VoxelsOnCartesianGrid
        Coord3DI image_size = {31, 111, 111};
        Coord3DF voxel_size = {3.375, 3, 3};
        IndexRange3D index_range(0, image_size.z() - 1,
                               -(image_size.y() / 2), -(image_size.y() / 2) + image_size.y() - 1,
                               -(image_size.x() / 2), -(image_size.x() / 2) + image_size.x() - 1);
        Coord3DF offset = {0.f, 0.f, 0.f};

        shared_ptr<Voxels3DF> im_sptr(new Voxels3DF(
            index_range,
			offset,
			voxel_size));
		im_sptr->fill(0.0);
        STIRImageData stir_im(im_sptr);

        // Test crop
        Coord3DI new_size = {3,2,5};
        Coord3DF zooms = {1.f,1.f,1.f};
        stir_im.zoom_image(zooms, offset, new_size, stir::ZoomOptions::preserve_sum);

        if (stir_im.dimensions()["z"] != new_size.at(1) ||
                stir_im.dimensions()["y"] != new_size.at(2) ||
                stir_im.dimensions()["x"] != new_size.at(3))
            throw std::runtime_error("STIRImageData::zoom_image failed");

#ifdef STIR_WITH_NiftyPET_PROJECTOR
        std::cout << "\nTesting NiftyPET projection...\n";
        // Load mMR sinogram
        const std::string f_mMR_template = f_template = path + "mMR_template_span11.hs";
        PETAcquisitionDataInFile mMR_template(f_mMR_template.c_str());
        std::shared_ptr<PETAcquisitionDataInMemory> acq_data_mMR_sptr(
                    new PETAcquisitionDataInMemory(mMR_template));
        acq_data_mMR_sptr->fill(1.f);

        // Create mMR image
        BasicCoordinate<3, int> min_image_indices(make_coordinate(0,  -160, -160));
        BasicCoordinate<3, int> max_image_indices(make_coordinate(126, 159,  159));
        IndexRange<3> range = IndexRange<3>(min_image_indices,max_image_indices);
        std::shared_ptr<STIRImageData> im_mMR_sptr(
                    new STIRImageData(Voxels3DF(
                                          acq_data_mMR_sptr->get_exam_info_sptr()->create_shared_clone(),
                                          range,
                                          CartesianCoordinate3D<float>(0.f,0.f,0.f),
                                          CartesianCoordinate3D<float>(2.03125f, 2.08626f, 2.08626f))));

        im_mMR_sptr->fill(1.f);
        PETAcquisitionModelUsingNiftyPET acq_model;
        stir::Verbosity::set(0);
        std::cout << "\nSetting up NiftyPET acquisition model...\n";
        acq_model.set_up(acq_data_mMR_sptr, im_mMR_sptr);
        std::cout << "\nForward projecting with NiftyPET acquisition model...\n";

        std::shared_ptr<PETAcquisitionData> prj_sptr = acq_model.forward(*im_mMR_sptr);
        std::cout << "\nBack projecting with NiftyPET acquisition model...\n";
        im_mMR_sptr = acq_model.backward(*prj_sptr);
        std::cout << "\nNiftyPET test succeeded.\n";
#endif

		std::cout << "done with test4.cpp...\n";
		return 0;
	}
	catch (...)
	{
		return 1;
	}
}
