#pragma once

#ifndef SIRF_ABSTRACT_PET_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_PET_IMAGE_DATA_TYPE

#include "sirf/common/ImageData.h"

namespace sirf {
	class PETImageData : public ImageData {
	public:
        /// Clone and return as unique pointer.
        std::unique_ptr<PETImageData> clone() const
        {
            return std::unique_ptr<PETImageData>(this->clone_impl());
        }
    protected:
        /// Clone helper function. Don't use.
        virtual PETImageData* clone_impl() const = 0;
	};
}

#endif