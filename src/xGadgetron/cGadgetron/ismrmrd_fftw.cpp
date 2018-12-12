/*
* copied from ismrmrd_fftw.h
*
*  Created on: Apr 1, 2013
*      Author: Michael S. Hansen

ISMRMRD SOFTWARE LICENSE JULY 2013

PERMISSION IS HEREBY GRANTED, FREE OF CHARGE, TO ANY PERSON OBTAINING
A COPY OF THIS SOFTWARE AND ASSOCIATED DOCUMENTATION FILES (THE
"SOFTWARE"), TO DEAL IN THE SOFTWARE WITHOUT RESTRICTION, INCLUDING
WITHOUT LIMITATION THE RIGHTS TO USE, COPY, MODIFY, MERGE, PUBLISH,
DISTRIBUTE, SUBLICENSE, AND/OR SELL COPIES OF THE SOFTWARE, AND TO
PERMIT PERSONS TO WHOM THE SOFTWARE IS FURNISHED TO DO SO, SUBJECT TO
THE FOLLOWING CONDITIONS:

THE ABOVE COPYRIGHT NOTICE, THIS PERMISSION NOTICE, AND THE LIMITATION
OF LIABILITY BELOW SHALL BE INCLUDED IN ALL COPIES OR REDISTRIBUTIONS
OF SUBSTANTIAL PORTIONS OF THE SOFTWARE.

SOFTWARE IS BEING DEVELOPED IN PART AT THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE, NATIONAL INSTITUTES OF HEALTH BY AN EMPLOYEE OF THE FEDERAL
GOVERNMENT IN THE COURSE OF HIS OFFICIAL DUTIES. PURSUANT TO TITLE 17,
SECTION 105 OF THE UNITED STATES CODE, THIS SOFTWARE IS NOT SUBJECT TO
COPYRIGHT PROTECTION AND IS IN THE PUBLIC DOMAIN. EXCEPT AS CONTAINED IN
THIS NOTICE, THE NAME OF THE AUTHORS, THE NATIONAL HEART, LUNG, AND BLOOD
INSTITUTE (NHLBI), OR THE NATIONAL INSTITUTES OF HEALTH (NIH) MAY NOT
BE USED TO ENDORSE OR PROMOTE PRODUCTS DERIVED FROM THIS SOFTWARE WITHOUT
SPECIFIC PRIOR WRITTEN PERMISSION FROM THE NHLBI OR THE NIH.THE SOFTWARE IS
PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include <fftw3.h>

#include "sirf/cGadgetron/ismrmrd_fftw.h"

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

