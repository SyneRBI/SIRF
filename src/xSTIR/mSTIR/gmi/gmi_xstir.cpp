/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

/*!
\file
\ingroup Matlab Interface Generators
\brief The executable for generating Matlab interface for cstir library.

\author Evgueni Ovtchinnikov
\author CCP PETMR
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
	int wp = 0);

int main()
{
	int status;
	string path_in;
	string path_out;
	string SIRF_path;
	SIRF_path = sirf::getenv("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		cout << "SIRF_PATH not defined, cannot find cstir library" << endl;
		return 1;
	}
	path_in = SIRF_path + "/src/xSTIR/cSTIR/include/";
	path_out = SIRF_path + "/src/xSTIR/mSTIR/";
	status = generate_matlab_interface
		("CSTIR", "cSTIR", 
			path_in, "sirf/STIR/cstir.h", 
			path_out, "mstir.h", "mstir.c", 1);
	if (status) {
		cout << "wrong input file format" << endl;
		return 1;
	}
}