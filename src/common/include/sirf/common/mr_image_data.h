#pragma once

#ifndef SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE

#include <ismrmrd/ismrmrd.h>

#include "image_data.h"

namespace sirf {
	class MRImageData : public ImageData<complex_float_t> {
	public:
	};
}

#endif