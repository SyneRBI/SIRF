/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

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

#include <iostream>
#include <string>

#include "sirf/Gadgetron/gadgetron_data_containers.h"
#include "sirf/STIR/stir_data_containers.h"
#include "sirf/Reg/NiftiImageData3D.h"
#include "sirf/Syn/utilities.h"

using namespace sirf;

int main(int argc, char* argv[])
{
	if (argc < 3) {
		std::cout << "usage: test_conv_img <filename> <engine>\n";
		return 1;
	}
	std::string filename(argv[1]);
	std::string eng_out(argv[2]);
	std::string eng_in;
	int len = filename.size();
	if (filename.substr(len - 3, 3) == std::string(".h5"))
		eng_in = "Gadgetron";
	else if (filename.substr(len - 3, 3) == std::string(".hv"))
		eng_in = "STIR";
	else if (filename.substr(len - 4, 4) == std::string(".nii"))
		eng_in = "Reg";
	else {
		std::cout << "unknown file format\n";
		return 1;
	}
	std::cout << "converting " << eng_in.c_str() << " image to "
		<< eng_out.c_str() << " image...\n";
	ImageDataWrap imw(filename, eng_in, true);
}
