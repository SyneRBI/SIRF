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
#include "sirf/common/iequals.h"
#include "sirf/common/utilities.h"

#include "object.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

int main()
{
	std::cout << "running test7.cpp...\n";
	try {
		std::string SIRF_data_path = examples_data_path("PET");
		if (SIRF_data_path.length() < 1) {
			std::cout << "cannot find data" << std::endl;
			return 1;
		}
		std::string path = append_path(SIRF_data_path, "mMR");
		std::cout << path << '\n';
		std::string f_listmode = append_path(path, "list.l.hdr");
		std::string f_template = append_path(path, "mMR_template_span11_small.hs");
		std::string f_mu_map = append_path(path, "mu_map.hv");
		
		STIRAcquisitionDataInMemory::set_as_template();

		//STIRAcquisitionDataInFile acq_data_template(f_template.c_str());
		CREATE_OBJECT(STIRAcquisitionData, STIRAcquisitionDataInFile,
			acq_data_template, templ_sptr, f_template.c_str());
		STIRListmodeData lm_data(f_listmode);
		//STIRImageData mu_map(f_mu_map);
		CREATE_OBJ(STIRImageData, mu_map, mu_map_sptr, f_mu_map);
		std::shared_ptr<STIRAcquisitionData> acq_sens_sptr;
		ListmodeToSinograms converter;

		std::shared_ptr<STIRAcquisitionData> sinograms_sptr;
		std::shared_ptr<STIRAcquisitionData> randoms_sptr;
		converter.sinograms_and_randoms_from_listmode
			(lm_data, 0, 10, acq_data_template, sinograms_sptr, randoms_sptr);
		
		std::cout << sinograms_sptr->norm() << '\n';
		std::cout << randoms_sptr->norm() << '\n';
		
		converter.acquisition_sensitivity_from_attenuation(templ_sptr, mu_map_sptr, acq_sens_sptr);
		std::cout << acq_sens_sptr->norm() << '\n';

		std::cout << "done with test7.cpp...\n";
		return 0;
	}
	catch (...)
	{
		return 1;
	}
}
