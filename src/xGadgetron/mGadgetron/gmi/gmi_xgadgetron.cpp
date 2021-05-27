/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
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
\ingroup MatlabInterface
\brief The executable for generating Matlab interface for cgadgetron library.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/

#include <stdio.h>

#include <iostream>

using namespace std;

#include "sirf/common/getenv.h"

// Matlab interface generator prototype - see src/common/mig.cpp
int generate_matlab_interface(
	const char* library,
	const char* prefix,
	const string& path_in,
	const string& chfile,
	const string& path_out,
	const string& mhfile,
	const string& mcfile,
	int wp = 0, bool to_cout = false);

int main(int argc, char *argv[])
{
	if (argc < 2) {
		cout << "Give output folder as input argument" << endl;
		return 1;
	}
	int status;
	string path_in;
	const string path_out = argv[1];
	string SIRF_path;
	if (argc >= 3)
		SIRF_path = argv[2];
	else
		SIRF_path = sirf::getenv("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		cout << "SIRF_PATH not defined, cannot find cgadgetron library" << endl;
		return 1;
	}
	path_in = SIRF_path + "/src/xGadgetron/cGadgetron/include/";
	status = generate_matlab_interface\
		("CGADGETRON", "cGT", 
			path_in, "sirf/Gadgetron/cgadgetron.h", 
			path_out, "/mgadgetron.h", "/mgadgetron.c");
	if (status) {
		cout << "wrong input file format" << endl;
		return 1;
	}
}
