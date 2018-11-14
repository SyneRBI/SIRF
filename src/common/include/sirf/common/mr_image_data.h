#pragma once

#ifndef SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE

#include <ismrmrd/ismrmrd.h>

#include "image_data.h"

namespace sirf {
	template<typename Iterator, typename Iterator_const>
	class MRImageData : public ImageData<complex_float_t> {
	public:
		virtual Iterator& begin() = 0;
		virtual Iterator_const& begin() const = 0;
		virtual Iterator& end() = 0;
		virtual Iterator_const& end() const = 0;
	};
	typedef std::iterator<std::forward_iterator_tag, NumRef> 
		ISMRMRDImageDataIterator;
	typedef std::iterator<std::forward_iterator_tag, NumRef const>
		ISMRMRDImageDataIterator_const;
	class ISMRMRDImageData : public MRImageData 
		< ISMRMRDImageDataIterator, ISMRMRDImageDataIterator_const > {
	};
	typedef std::iterator<std::forward_iterator_tag, std::complex<float> >
		ComplexImageDataIterator;
	typedef std::iterator<std::forward_iterator_tag, std::complex<float> const >
		ComplexImageDataIterator_const;
	class ComplexImageData : public MRImageData 
		< ComplexImageDataIterator, ComplexImageDataIterator_const > {
	};
}

#endif