#pragma once

#include "sirf_image_data.h"

namespace sirf {
	template <typename T>
	class PETImageData : public SIRFImageData<T>
	{
	public:
		virtual ~PETImageData() {}
	};
}