#pragma once

#ifndef SIRF_ABSTRACT_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_IMAGE_DATA_TYPE

#include "data_container.h"
#include "geometrical_info.h"

/*!
\ingroup SIRFImageDataClasses
\brief Abstract base class for SIRF image data.

*/
namespace sirf {
	class ImageData : public DataContainer
	{
	public:
		virtual ~ImageData() {}
		virtual Dimensions dimensions() const = 0; // to go to DataContainer eventually
		//virtual void get_data(void* data) const = 0;
		//virtual void set_data(const void* data) = 0;
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
		void copy(Iterator_const& src, Iterator& dst, Iterator& end)
		{
			for (; dst != end; ++dst, ++src)
				*dst = *src;
		}
        /// Get geometrical info
        std::shared_ptr<const VoxelisedGeometricalInfo3D > get_geom_info()
        {
            // If the geometrical info has not been created yet, do that
            if (!_geom_info_sptr)
                this->set_up_geom_info();
            return _geom_info_sptr;
        }
    protected:
        /// Populate the geometrical info metadata (from the image's own metadata)
        virtual void set_up_geom_info() = 0;
        std::shared_ptr<VoxelisedGeometricalInfo3D> _geom_info_sptr;
	};
}

#endif
