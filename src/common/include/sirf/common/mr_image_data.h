#pragma once

#ifndef SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE

#include <ismrmrd/ismrmrd.h>

#include "image_data.h"

namespace sirf {
	template<typename ImageDataIterator>
	class MRImageData : public ImageData<complex_float_t> {
	public:
		virtual ImageDataIterator& begin() = 0;
		virtual ImageDataIterator& begin() const = 0;
		virtual ImageDataIterator& end() = 0;
		virtual ImageDataIterator& end() const = 0;
	};
}

#endif