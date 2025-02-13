/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2023 Rutherford Appleton Laboratory STFC
Copyright 2023 University College London

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
\ingroup Common

\author Evgueni Ovtchinnikov
\author Kris Thielemans
*/
#include "sirf/common/utilities.h"
#include "sirf/common/version.h"
#include "sirf/common/getenv.h"
#include <sstream>

namespace sirf {
	char path_separator()
	{
		std::string filename = __FILE__;
		auto found = filename.find_last_of("/\\");
		if (found == std::string::npos)
			return '/';
		return filename[found];
	}

	std::string examples_data_path(const char* data_type)
	{
		std::string SIRF_data_path = sirf::getenv("SIRF_DATA_PATH");
		if (SIRF_data_path.length() > 0)
			return append_path(SIRF_data_path, "examples", data_type);
		std::string SIRF_install_path = sirf::getenv("SIRF_INSTALL_PATH");
		if (SIRF_install_path.length() > 0) {
			std::stringstream sirf_version;
			sirf_version << "SIRF-" << SIRF_VERSION_MAJOR << '.' << SIRF_VERSION_MINOR;
			return append_path(SIRF_install_path, "share", sirf_version.str(), "data", "examples", data_type);
		}
		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() > 0)
			return append_path(SIRF_path, "data", "examples", data_type);
		return "";
	}
}
