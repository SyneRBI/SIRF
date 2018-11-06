#pragma once

#include "image_data.h"

namespace sirf {
	//template <typename T>
	class PETImageData : public ImageData<float>
	{
	public:
		virtual void get_data(float* data) const = 0;
		virtual void set_data(const float* data) = 0;
	};
}