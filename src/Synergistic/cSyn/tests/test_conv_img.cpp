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
	if (argc < 4) {
		std::cout << "usage: test_conv_img <filename> <engine_in> <engine_out>\n";
		return 1;
	}
	std::string filename(argv[1]);
	std::string eng_in(argv[2]);
	std::string eng_out(argv[3]);
	std::cout << "creating " << eng_in.c_str() << " image\n";
	ImageDataWrap imw(filename, eng_in, true);
	std::cout << "converting " << eng_in.c_str() << " image to "
		<< eng_out.c_str() << " image...\n";
	const ImageData& im_in = imw.data();
	if (eng_out == std::string("Reg")) {
		NiftiImageData3D<float> im_out(im_in);
		if (im_out == im_in) {
			std::cout << "images are identical\n";
			return 0;
		}
		else {
			std::cout << "images are not identical\n";
			return 1;
		}
	}
	else
		std::cout << "engine " << eng_out << " not supported yet\n";
	return 1;
}