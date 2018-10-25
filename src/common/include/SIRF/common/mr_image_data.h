#pragma once

#include "sirf_image_data.h"

namespace sirf {
	template <typename T>
	class MRImageData : public SIRFImageData<T>
	{
	public:
		virtual ~MRImageData() {}
	};
}