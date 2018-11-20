#pragma once

#ifndef SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE

#include <ismrmrd/ismrmrd.h>

#include "image_data.h"
#include "num_ref.h"

namespace sirf {
	//template<typename Iterator, typename Iterator_const>
	class MRImageData : public ImageData<complex_float_t> {
	public:
		class Iterator {
		public:
			virtual ~Iterator() {}
			virtual Iterator& operator++() = 0;
			virtual aNumRef& operator*() = 0;
			virtual bool operator==(const Iterator&) const = 0;
			virtual bool operator!=(const Iterator&) const = 0;
		};
		class Iterator_const {
		public:
			virtual ~Iterator_const() {}
			virtual Iterator_const& operator++() = 0;
			virtual const aNumRef& operator*() const = 0;
			virtual bool operator==(const Iterator_const&) const = 0;
			virtual bool operator!=(const Iterator_const&) const = 0;
		};
		virtual Iterator& begin() = 0;
		virtual Iterator_const& begin() const = 0;
		virtual Iterator& end() = 0;
		virtual Iterator_const& end() const = 0;
	};
}

#endif