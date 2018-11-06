#pragma once

#ifndef SIRF_ABSTRACT_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_IMAGE_DATA_TYPE

#include "data_container.h"

/*!
\ingroup SIRFImageDataClasses
\brief Abstract base class for SIRF image data.

*/
namespace sirf {
	template<typename T>
	class ImageData : public aDataContainer<T>
	{
	public:
		virtual ~ImageData() {}
		virtual void get_data(T* data) const = 0;
		virtual void set_data(const T* data) = 0;
	};
}

#endif