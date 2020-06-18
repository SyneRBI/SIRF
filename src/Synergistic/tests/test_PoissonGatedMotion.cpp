/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 University College London

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
\ingroup STIR Tests

\author Richard Brown
\author SyneRBI
*/

// STIR stuff
#include "sirf/STIR/stir_x.h"
#include "sirf/common/getenv.h"

// Reg stuff
#include "sirf/Reg/NiftiImageData3DDisplacement.h"

using namespace sirf;
using namespace stir;

int test4()
{
	try {

		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() < 1) {
			std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
			return 1;
		}

        // Open acquisition data
		std::string path = SIRF_path + "/data/examples/PET/";
		std::string f_acq_data = path + "Utahscat600k_ca_seg4.hs";
        auto acq_data_sptr = std::make_shared<PETAcquisitionDataInMemory>(f_acq_data.c_str());

        // Recon image
        auto im_sptr = std::make_shared<STIRImageData>(*acq_data_sptr);
        im_sptr->fill(0.f);

        // Create blank displacement field and convert to STIR format
        NiftiImageData3DDisplacement<float> disp;
        disp.create_from_3D_image(*im_sptr);
        const std::string temp_fname = "temp_disp.nii";
        disp.write(temp_fname);
        auto disp_sptr = MAKE_SHARED<STIRDisp3DF>(temp_fname, 3 /*bspline_order*/);

        // create additive term
		auto a_sptr = acq_data_sptr->new_acquisition_data();
		a_sptr->fill(0.05f);
		// create background term
		auto b_sptr = acq_data_sptr->new_acquisition_data();
		b_sptr->fill(0.1f);

        // Create acquisition model
        auto matrix_sptr = MAKE_SHARED<RayTracingMatrix>();
        matrix_sptr->set_num_tangential_LORs(2);
        auto am_sptr = MAKE_SHARED<AcqModUsingMatrix3DF>();
        am_sptr->set_matrix(matrix_sptr);
        am_sptr->set_additive_term(a_sptr);
        am_sptr->set_background_term(b_sptr);
        am_sptr->set_up(acq_data_sptr, im_sptr);

        // Create poisson loglikelihood
        PoissonLogLhLinModMeanGatedProjDataWMotion3DF obj_fn;

        // Add motion gates
        const unsigned num_motion_gates = 4;
        for (unsigned i=0; i<num_motion_gates; ++i)
            obj_fn.add_gate(acq_data_sptr, am_sptr, disp_sptr);

        obj_fn.set_up(im_sptr->data_sptr());

        // Recon
        OSMAPOSLReconstruction<Image3DF> recon;
        recon.set_up(im_sptr->data_sptr());
        recon.set_num_subiterations(2);
        recon.reconstruct(im_sptr->data_sptr());

        std::cout << "\nSuccess!\n";

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
