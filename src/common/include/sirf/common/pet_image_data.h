#pragma once

#ifndef SIRF_ABSTRACT_PET_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_PET_IMAGE_DATA_TYPE

#include "image_data.h"

namespace sirf {
	template<typename Iterator, typename Iterator_const>
	class PETImageData : public ImageData<float> {
	public:
		virtual Iterator& begin() = 0;
		virtual Iterator_const& begin() const = 0;
		virtual Iterator& end() = 0;
		virtual Iterator_const& end() const = 0;
	};
}

#endif