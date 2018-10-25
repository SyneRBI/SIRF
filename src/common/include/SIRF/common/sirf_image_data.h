#pragma once

#include "data_container.h"

/*
\ingroup SIRF Image Data classes
\brief Abstract base class for SIRF image geometry and data.

*/
namespace sirf {
	template<typename T>
	class SIRFImageData : public aDataContainer<T>
	{
	public:
		virtual ~SIRFImageData() {}
	};
}