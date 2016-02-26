#ifndef CGADGETRON_ISMRMRD_FFTW
#define CGADGETRON_ISMRMRD_FFTW
/*
* ismrmrd_fftw.h
*
*  Created on: Apr 1, 2013
*      Author: Michael S. Hansen
*/

namespace ISMRMRD {
	template<typename TI, typename TO> 
	void 
		circshift(TO *out, const TI *in, int xdim, int ydim, int xshift, int yshift)
	{
		for (int i = 0; i < ydim; i++) {
			int ii = (i + yshift) % ydim;
			for (int j = 0; j < xdim; j++) {
				int jj = (j + xshift) % xdim;
				out[ii * xdim + jj] = in[i * xdim + j];
			}
		}
	}
	int fft2c(NDArray<complex_float_t> &a);
	int ifft2c(NDArray<complex_float_t> &a);

};

#endif