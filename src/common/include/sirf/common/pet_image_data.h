#pragma once

#ifndef SIRF_ABSTRACT_PET_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_PET_IMAGE_DATA_TYPE

#include "image_data.h"
#include "num_ref.h"

namespace sirf {
	//template<typename Iterator, typename Iterator_const>
	class PETImageData : public ImageData<float> {
	public:
		class Iter {
		public:
			virtual ~Iter() {}
			virtual Iter& operator++() = 0;
			virtual Iter& operator++(int) = 0;
			virtual aNumRef& operator*() = 0;
			virtual bool operator==(const Iter&) const = 0;
			virtual bool operator!=(const Iter&) const = 0;
		};
		class Iter_const {
		public:
			virtual ~Iter_const() {}
			virtual Iter_const& operator++() = 0;
			virtual Iter_const& operator++(int) = 0;
			virtual const aNumRef& operator*() const = 0;
			virtual bool operator==(const Iter_const&) const = 0;
			virtual bool operator!=(const Iter_const&) const = 0;
		};
		virtual Iter& begin_new() = 0;
		virtual Iter_const& begin_new() const = 0;
		virtual Iter& end_new() = 0;
		virtual Iter_const& end_new() const = 0;
		//virtual Iterator& begin() = 0;
		//virtual Iterator_const& begin() const = 0;
		//virtual Iterator& end() = 0;
		//virtual Iterator_const& end() const = 0;
	};
}

#endif