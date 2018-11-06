#pragma once

#include <ismrmrd/ismrmrd.h>

#include "image_data.h"

namespace sirf {
	//template <typename T>
	class MRImageData : public ImageData<complex_float_t> {
	public:
		virtual void get_data(complex_float_t* data) const = 0;
		virtual void set_data(const complex_float_t* data) = 0;
	};
}