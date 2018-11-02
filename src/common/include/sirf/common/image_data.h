#pragma once

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
	};
}