#pragma once

#ifndef SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_MR_IMAGE_DATA_TYPE

#include <ismrmrd/ismrmrd.h>

#include "sirf/common/ImageData.h"

namespace sirf {
	class MRImageData : public ImageData {
	public:
        /// Clone and return as unique pointer.
        std::unique_ptr<MRImageData> clone() const
        {
            return std::unique_ptr<MRImageData>(this->clone_impl());
        }
    protected:
        /// Clone helper function. Don't use.
        virtual MRImageData* clone_impl() const = 0;
	};
}

#endif