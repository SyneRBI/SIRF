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

int test1();
int test2();

int main()
{
	int failed = 0;
	failed += test1();
	failed += test2();
	std::cout << failed << " tests failed\n";
	return failed;
}
