/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2019 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2019 - 2020 University College London

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
\ingroup Gadgetron Extensions
\brief MR related C++ tests

\author Johannes Mayer
\author SyneRBI
*/

#include <iostream>
#include <cstdlib>

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/Gadgetron/gadgetron_x.h"
#include "sirf/Gadgetron/chain_lib.h"

#include "mrtest_auxiliary_funs.h"

using namespace sirf;

bool test_get_kspace_order(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av_slice;
        av_slice.read(fname_input);
        av_slice.sort();

        auto kspace_sorting_slice = av_slice.get_kspace_order();

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_get_subset(const std::string& fname_input)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        sirf::AcquisitionsVector av;
        av.read(fname_input);

        std::vector<int> subset_idx;
        for(int i=0; i<av.number()/10; ++i)
            subset_idx.push_back(i);

        sirf::AcquisitionsVector subset;
        av.get_subset(subset, subset_idx);

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_CoilSensitivitiesVector_calculate(const MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        CoilSensitivitiesVector csv;
        csv.set_csm_smoothness(50);
        csv.calculate(av);

        std::cout << "We have " << csv.items() << " coilmaps" << std::endl;

        for(int i=0; i<csv.items(); ++i)
        {
            gadgetron::shared_ptr<ImageWrap> sptr_iw = csv.sptr_image_wrap(i);

            std::stringstream fname_out;
            fname_out << "output_" << __FUNCTION__ << "_" << i;

            sirf::write_cfimage_to_raw(fname_out.str(), *sptr_iw);
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_CoilSensitivitiesVector_get_csm_as_cfimage(const MRAcquisitionData& av)
{
    try
    {
        std::cout << "Running test " << __FUNCTION__ << std::endl;

        CoilSensitivitiesVector csv;
        csv.calculate(av);

        std::cout << "We have " << csv.items() << " coilmaps" << std::endl;

        for(int i=0; i<csv.items(); ++i)
        {
            CFImage img = csv.get_csm_as_cfimage(i);

            std::stringstream fname_out;
            fname_out << "output_" << __FUNCTION__ << "_" << i;

            sirf::write_cfimage_to_raw(fname_out.str(), img);
        }

        return true;

    }
    catch( std::runtime_error const &e)
    {
        std::cout << "Exception caught " <<__FUNCTION__ <<" .!" <<std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
}

bool test_acq_mod_norm(gadgetron::shared_ptr<MRAcquisitionData> sptr_ad)
{
    MRAcquisitionData& ad = *sptr_ad;
    int acq_dim[10];
    ad.get_acquisitions_dimensions((size_t)acq_dim);
    std::cout << "acquisitions dimensions: "
              << acq_dim[0] << " by "
              << acq_dim[1] << " by "
              << acq_dim[2] << '\n';

    gadgetron::shared_ptr<CoilSensitivitiesVector> sptr_csv
        (new CoilSensitivitiesVector);
    CoilSensitivitiesVector& csv = *sptr_csv;
    csv.calculate(ad);

    gadgetron::shared_ptr<GadgetronImageData> sptr_id;
    gadgetron::shared_ptr<GadgetronImageData> sptr;
    if (ad.undersampled()) {
        std::cout << "found undersampled data\n";
        SimpleGRAPPAReconstructionProcessor recon;
        std::cout << "reconstructing...\n";
        recon.process(ad);
        sptr_id = recon.get_output();
        sptr_id = sptr_id->clone("GADGETRON_DataRole", "image");
    }
    else {
        std::cout << "found fully sampled data\n";
        SimpleReconstructionProcessor recon;
        std::cout << "reconstructing...\n";
        recon.process(ad);
        sptr_id = recon.get_output();
    }

    GadgetronImageData& image = *sptr_id;
    int img_dim[10];
    int ni = image.number();
    image.get_image_dimensions(0, img_dim);
    std::cout << ni << " images reconstructed\n";
    std::cout << "image dimensions: "
              << img_dim[0] << " by "
              << img_dim[1] << " by "
              << img_dim[2] << " by "
              << img_dim[3] << '\n';

    MRAcquisitionModel am;
    am.set_up(sptr_ad, sptr_id);
    am.setCSMs(sptr_csv);
    std::cout << "forward projectiong...\n";
    gadgetron::shared_ptr<MRAcquisitionData> sptr_sd = am.fwd(image);
    MRAcquisitionData& sd = *sptr_sd;
    sd.get_acquisitions_dimensions((size_t)acq_dim);
    std::cout << "simulated acquisitions dimensions: "
              << acq_dim[0] << " by "
              << acq_dim[1] << " by "
              << acq_dim[2] << '\n';

    float im_norm = image.norm();
    float sd_norm = sd.norm();
    float am_norm = am.norm();
    std::cout << "\nchecking the acquisition model norm:\n";
    std::cout << "acquisition model norm: |A| = " << am_norm << '\n';
    std::cout << "image data x norm: |x| = " << im_norm << '\n';
    std::cout << "simulated acquisition data norm: |A(x)| = " << sd_norm << '\n';
    std::cout << "checking that |A(x)| <= |A||x|: ";
    float bound = am_norm*im_norm;
    bool ok = (sd_norm <= bound);
    if (ok)
        std::cout << sd_norm << " <= " << bound << " ok!\n";
    else
        std::cout << sd_norm << " > " << bound << " failure!\n";

    return ok;
}

#include <ismrmrd/ismrmrd.h>

void edit_acq(const char* in, const char* out)
{
	std::string SIRF_PATH = getenv("SIRF_PATH");
	std::string path_in = SIRF_PATH + in;
	std::string path_out = out;
	gadgetron::shared_ptr<MRAcquisitionData> sptr_ad(new AcquisitionsVector);
	AcquisitionsVector& av = (AcquisitionsVector&)*sptr_ad;
	av.read(path_in);
	int na = av.number();
	std::cout << na << " acquisitions read from " << path_in << '\n';
	int acq_dim[10];
	av.get_acquisitions_dimensions((size_t)acq_dim);
	std::cout << "acquisitions dimensions: "
		<< acq_dim[0] << " by "
		<< acq_dim[1] << " by "
		<< acq_dim[2] << '\n';
	ISMRMRD::Acquisition acq;
	for (int i = 0; i < na; i++) {
		av.get_acquisition(i, acq);
		float* read_dir = acq.read_dir();
		float* phase_dir = acq.phase_dir();
		float* slice_dir = acq.slice_dir();
		if (i == 0) {
			std::cout << "in:\n";
			std::cout << read_dir[0] << ',' << read_dir[1] << ',' << read_dir[2] << '\n';
			std::cout << phase_dir[0] << ',' << phase_dir[1] << ',' << phase_dir[2] << '\n';
			std::cout << slice_dir[0] << ',' << slice_dir[1] << ',' << slice_dir[2] << '\n';
		}
		read_dir[0] = 1.0f;
		read_dir[1] = 0.0f;
		read_dir[2] = 0.0f;
		phase_dir[0] = 0.0f;
		phase_dir[1] = 1.0f;
		phase_dir[2] = 0.0f;
		slice_dir[0] = 0.0f;
		slice_dir[1] = 0.0f;
		slice_dir[2] = 1.0f;
		if (i == 0) {
			std::cout << "out:\n";
			std::cout << read_dir[0] << ',' << read_dir[1] << ',' << read_dir[2] << '\n';
			std::cout << phase_dir[0] << ',' << phase_dir[1] << ',' << phase_dir[2] << '\n';
			std::cout << slice_dir[0] << ',' << slice_dir[1] << ',' << slice_dir[2] << '\n';
		}
		av.set_acquisition(i, acq);
	}
	std::cout << "saving edited acquisition data to " << out << "...";
	av.write(out);
	std::cout << "done\n";
}

int main ( int argc, char* argv[])
{

	try{
		if (argc == 3) {
			edit_acq(argv[1], argv[2]);
			return 0;
		}

        std::string SIRF_PATH;
        if (argc==1)
            SIRF_PATH = getenv("SIRF_PATH");
        else
            SIRF_PATH = argv[1];

        std::string data_path = SIRF_PATH + "/data/examples/MR/simulated_MR_2D_cartesian_Grappa2.h5";

//        test_get_kspace_order(data_path);
//        test_get_subset(data_path);

//        sirf::AcquisitionsVector av;
        gadgetron::shared_ptr<MRAcquisitionData> sptr_ad(new AcquisitionsVector);
        AcquisitionsVector& av = (AcquisitionsVector&)*sptr_ad;
        av.read(data_path);

        sirf::preprocess_acquisition_data(av);

        test_CoilSensitivitiesVector_calculate(av);
        test_CoilSensitivitiesVector_get_csm_as_cfimage(av);
        test_acq_mod_norm(sptr_ad);
        return 0;
	}
    catch(const std::exception &error) {
        std::cerr << "\nHere's the error:\n\t" << error.what() << "\n\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

