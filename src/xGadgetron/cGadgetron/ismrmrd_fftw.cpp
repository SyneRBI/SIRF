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

#include <complex>
#include <iostream>
#include <mutex>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include <fftw3.h>

#include "sirf/Gadgetron/ismrmrd_fftw.h"

typedef complex_float_t ComplexType;

#define USE_OMP

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

	fftwf_plan fftw_plan_dft_3d_(int n0, int n1, int n2, ComplexType* in, ComplexType * out, int sign, unsigned int flags) {
		return fftwf_plan_dft_3d(n0, n1, n2, (fftwf_complex*)in, (fftwf_complex*)out, sign, flags);
	}

	void fftw_execute_dft_(fftwf_plan_s * ptr, ComplexType* in, ComplexType* out) {
		fftwf_execute_dft(ptr, (fftwf_complex*)in, (fftwf_complex*)out);
	}

	void fftw_destroy_plan_(fftwf_plan p) {
		fftwf_destroy_plan(p);
	}

	void fftshiftPivot3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty, size_t pivotz)
	{

		if (a == NULL)
			throw std::runtime_error("fftshiftPivot3D: void ptr provided");

		long long tt;

#pragma omp parallel private(tt) shared(a, x, y, z, n, pivotx, pivoty, pivotz) if (n>16)
		{
			//hoNDArray< ComplexType > aTmp(x*y*z);
			ComplexType* tmp =
				(ComplexType*)fftwf_malloc(sizeof(ComplexType)*x*y*z);

#pragma omp for
			for (tt = 0; tt < (long long)n; tt++)
			{
				size_t ay, ry, az, rz;

				for (az = pivotz; az < z; az++)
				{
					rz = az - pivotz;

					const ComplexType* ac = a + tt*x*y*z + az*x*y;
					ComplexType* rc = tmp + rz*x*y;
					//ComplexType* rc = aTmp.begin() + rz*x*y;

					for (ay = pivoty; ay < y; ay++)
					{
						ry = ay - pivoty;
						memcpy(rc + ry*x, ac + ay*x + pivotx, sizeof(ComplexType)*(x - pivotx));
						memcpy(rc + ry*x + x - pivotx, ac + ay*x, sizeof(ComplexType)*pivotx);
					}

					for (ay = 0; ay < pivoty; ay++)
					{
						ry = ay + y - pivoty;
						memcpy(rc + ry*x, ac + ay*x + pivotx, sizeof(ComplexType)*(x - pivotx));
						memcpy(rc + ry*x + x - pivotx, ac + ay*x, sizeof(ComplexType)*pivotx);
					}
				}

				for (az = 0; az < pivotz; az++)
				{
					rz = az + z - pivotz;

					const ComplexType* ac = a + tt*x*y*z + az*x*y;
					ComplexType* rc = tmp + rz*x*y;
					//ComplexType* rc = aTmp.begin() + rz*x*y;

					for (ay = pivoty; ay < y; ay++)
					{
						ry = ay - pivoty;
						memcpy(rc + ry*x, ac + ay*x + pivotx, sizeof(ComplexType)*(x - pivotx));
						memcpy(rc + ry*x + x - pivotx, ac + ay*x, sizeof(ComplexType)*pivotx);
					}

					for (ay = 0; ay < pivoty; ay++)
					{
						ry = ay + y - pivoty;
						memcpy(rc + ry*x, ac + ay*x + pivotx, sizeof(ComplexType)*(x - pivotx));
						memcpy(rc + ry*x + x - pivotx, ac + ay*x, sizeof(ComplexType)*pivotx);
					}
				}

				memcpy(a + tt*x*y*z, tmp, sizeof(ComplexType)*x*y*z);
				//memcpy(a + tt*x*y*z, aTmp.begin(), sizeof(ComplexType)*x*y*z);
			}
		}
	}

	inline size_t fftshiftPivot(size_t x)
	{
		return (size_t)(ceil(x*0.5));
	}

	inline size_t ifftshiftPivot(size_t x)
	{
		return (size_t)(floor(x*0.5));
	}

	inline void fftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n)
	{
		if (a == NULL)
			throw std::runtime_error("fftshift3D: void ptr provided");
		size_t pivotx = fftshiftPivot(x);
		size_t pivoty = fftshiftPivot(y);
		size_t pivotz = fftshiftPivot(z);
		fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz);
	}

	inline void ifftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n)
	{
		if (a == NULL) throw std::runtime_error("hoNDFFT::ifftshift3D: void ptr provided");

		size_t pivotx = ifftshiftPivot(x);
		size_t pivoty = ifftshiftPivot(y);
		size_t pivotz = ifftshiftPivot(z);

		fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz);
	}

	inline void fftshift3D(NDArray< ComplexType >& a)
	{
		const size_t* dims = a.getDims();
		size_t n = a.getNumberOfElements() / (dims[0] * dims[1] * dims[2]);
		return fftshift3D(a.begin(), dims[0], dims[1], dims[2], n);
	}

	inline void ifftshift3D(NDArray< ComplexType >& a)
	{
		const size_t* dims = a.getDims();
		size_t n = a.getNumberOfElements() / (dims[0] * dims[1] * dims[2]);
		return ifftshift3D(a.begin(), dims[0], dims[1], dims[2], n);
	}

	inline int get_num_threads_fft3(size_t n0, size_t n1, size_t n2, size_t num)
	{
		int num_of_max_threads_;
#ifdef USE_OMP
		num_of_max_threads_ = 8; //omp_get_num_procs();
#else
		num_of_max_threads_ = 1;
#endif // USE_OMP
		if (num_of_max_threads_ == 1)
			return 1;

		if (num >= num_of_max_threads_) {
			return num_of_max_threads_;
		}

		return 1;
	}

	void fft3(NDArray< ComplexType >& a, NDArray< ComplexType >& r, bool forward)
	{
		r = a;

		const size_t* dims = a.getDims();
		int n2 = (int)dims[0];
		int n1 = (int)dims[1];
		int n0 = (int)dims[2];

		float fftRatio = float(1.0 / std::sqrt(float(n0*n1*n2)));

		int num = (int)(a.getNumberOfElements() / (n0*n1*n2));
		int num_thr = get_num_threads_fft3(n0, n1, n2, num);

		long long n;

		std::mutex mutex_;
		fftwf_plan p;
		//typename fftw_types<T>::plan * p;

		{
			std::lock_guard<std::mutex> guard(mutex_);
			p = fftw_plan_dft_3d_(n0, n1, n2,
				a.getDataPtr(),
				r.getDataPtr(),
				forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

		}

#pragma omp parallel for private(n) shared(num, p, a, n0, n1, n2, r) if (num_thr > 1) num_threads(num_thr)
		for (n = 0; n < num; n++)
		{
			fftw_execute_dft_(p, a.begin() + n*n0*n1*n2,
				r.begin() + n*n0*n1*n2);
		}

		{
			std::lock_guard<std::mutex> guard(mutex_);
			fftw_destroy_plan_(p);
		}

		for (size_t n = 0; n < a.getNumberOfElements(); n++) {
			r.getDataPtr()[n] *= fftRatio;
		}
		//	r *= fftRatio;

	}

	void fft3(NDArray< ComplexType >& a, bool forward)
	{
		NDArray< ComplexType > res(a);
		fft3(res, a, forward);
	}

	inline void fft3(NDArray< ComplexType >& a)
	{
		fft3(a, true);
	}

	inline void ifft3(NDArray< ComplexType >& a)
	{
		fft3(a, false);
	}

	void fft3c(NDArray< ComplexType >& a)
	{
		ifftshift3D(a);
		fft3(a);
		fftshift3D(a);
	}

	void ifft3c(NDArray< ComplexType >& a)
	{
		ifftshift3D(a);
		ifft3(a);
		fftshift3D(a);
	}

}