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
\ingroup Gadgetron Image Wrapper
\brief Iplementation file for a wrapper class for ISMRMRD::Image.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#include "gadgetron_image_wrap.h"

void
ImageWrap::get_cmplx_data(float* re, float* im) const
{
	int dim[4];
	size_t n = get_dim(dim);
	if (type_ == ISMRMRD::ISMRMRD_CXFLOAT) {
		const CFImage& img = *(const CFImage*)ptr_;
		const complex_float_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++) {
			complex_float_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
	}
	else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE) {
		const CDImage& img = *(const CDImage*)ptr_;
		const complex_double_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++) {
			complex_double_t z = ptr[i];
			re[i] = std::real(z);
			im[i] = std::imag(z);
		}
	}
	else {
		get_data(re);
		for (size_t i = 0; i < n; i++)
			im[i] = 0;
	}
}

void
ImageWrap::set_cmplx_data(const float* re, const float* im) const
{
	int dim[4];
	size_t n = get_dim(dim);
	if (type_ == ISMRMRD::ISMRMRD_CXFLOAT) {
		CFImage& img = *(CFImage*)ptr_;
		complex_float_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++)
			ptr[i] = std::complex<float>((float)re[i], (float)im[i]);
	}
	else if (type_ == ISMRMRD::ISMRMRD_CXDOUBLE) {
		CDImage& img = *(CDImage*)ptr_;
		complex_double_t* ptr = img.getDataPtr();
		for (size_t i = 0; i < n; i++)
			ptr[i] = std::complex<double>(re[i], im[i]);
	}
}
