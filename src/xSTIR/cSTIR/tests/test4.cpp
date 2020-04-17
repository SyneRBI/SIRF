#include <iostream>
#include <cstdlib>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"

#include "sirf/STIR/stir_x.h"
#include "sirf/common/getenv.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

int test4()
{
	try {

		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

		std::string path = SIRF_path + "/data/examples/PET/mMR/";

		string f_listmode = path + "list.l.hdr";
		string f_template = path + "mMR_template_span11_small.hs";

		// Listmode to sinograms
		ListmodeToSinograms converter;
		converter.set_input(f_listmode);
		converter.set_output("proj_data");
		converter.set_template(f_template);
		converter.set_time_interval(0, 10);
		converter.set_up();
		converter.estimate_randoms();

        // Check count rates - for the particular dataset,
        // we know that 73036 is exceeded at 22s. You can
        // see this with STIR's list_lm_countrates
        const float prompt_rate_threshold = 73036.f;
        const float known_time = 22.f;
        const float time_at_which_prompt_rate_exceeds_threshold =
                converter.get_time_at_which_prompt_rate_exceeds_threshold(prompt_rate_threshold);
        if (std::abs(time_at_which_prompt_rate_exceeds_threshold-known_time) > 1e-4f)
            throw std::runtime_error("ListmodeToSinograms::get_time_at_which_prompt_rate_exceeds_threshold failed");

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

		return 0;
	}
	catch (...)
	{
		return 1;
	}
}

//int test5();

int main()
{
	return test4();
	//return test5();
}