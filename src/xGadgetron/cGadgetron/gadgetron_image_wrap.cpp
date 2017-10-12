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
