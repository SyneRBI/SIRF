/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
Copyright 2015 - 2017 University College London.

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

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include <fftw3.h>

#include "ismrmrd_fftw.h"

namespace ISMRMRD {

#define fftshift(out, in, x, y) circshift(out, in, x, y, (x/2), (y/2))

	int fft2c(NDArray<complex_float_t> &a, bool forward)
	{
		if (a.getNDim() < 2) {
			std::cout << "fft2c Error: input array must have at least two dimensions"
				<< std::endl;
			return -1;
		}

		size_t elements = a.getDims()[0] * a.getDims()[1];
		size_t ffts = a.getNumberOfElements() / elements;

		//Array for transformation
		fftwf_complex* tmp =
			(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*a.getNumberOfElements());

		if (!tmp) {
			std::cout << "Error allocating temporary storage for FFTW" << std::endl;
			return -1;
		}

		for (size_t f = 0; f < ffts; f++) {

			fftshift(reinterpret_cast<std::complex<float>*>(tmp),
				&a(0, 0, f), a.getDims()[0], a.getDims()[1]);

			//Create the FFTW plan
			fftwf_plan p;
			if (forward) {
				p = fftwf_plan_dft_2d
					(a.getDims()[1], a.getDims()[0], tmp, tmp,
					FFTW_FORWARD, FFTW_ESTIMATE);
			}
			else {
				p = fftwf_plan_dft_2d
					(a.getDims()[1], a.getDims()[0], tmp, tmp,
					FFTW_BACKWARD, FFTW_ESTIMATE);
			}
			fftwf_execute(p);

			fftshift(&a(0, 0, f), reinterpret_cast<std::complex<float>*>(tmp),
				a.getDims()[0], a.getDims()[1]);

			//Clean up.
			fftwf_destroy_plan(p);
		}

		std::complex<float> scale(std::sqrt(1.0f*elements), 0.0);
		for (size_t n = 0; n < a.getNumberOfElements(); n++) {
			a.getDataPtr()[n] /= scale;
		}
		fftwf_free(tmp);
		return 0;
	}

	int fft2c(NDArray<complex_float_t> &a) 
	{
		return fft2c(a, true);
	}

	int ifft2c(NDArray<complex_float_t> &a) 
	{
		return fft2c(a, false);
	}

};

