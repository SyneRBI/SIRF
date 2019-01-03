#pragma once

#ifndef SIRF_ABSTRACT_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_IMAGE_DATA_TYPE

#include "sirf/common/ANumRef.h"
#include "sirf/common/DataContainer.h"
#include "sirf/common/DataContainer.h"
#include "sirf/common/ANumRef.h"
#include "sirf/common/GeometricalInfo.h"

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
                        virtual ANumRef& operator*() = 0;
			virtual bool operator==(const Iterator&) const = 0;
			virtual bool operator!=(const Iterator&) const = 0;
		};
		class Iterator_const {
		public:
			virtual ~Iterator_const() {}
			virtual Iterator_const& operator++() = 0;
                        virtual const ANumRef& operator*() const = 0;
			virtual bool operator==(const Iterator_const&) const = 0;
			virtual bool operator!=(const Iterator_const&) const = 0;
		};
		virtual Iterator& begin() = 0;
		virtual Iterator_const& begin() const = 0;
		virtual Iterator& end() = 0;
		virtual Iterator_const& end() const = 0;
		virtual bool ordered() const
		{
			return true;
		}
		void copy(Iterator_const& src, Iterator& dst, Iterator& end) const
		{
			for (; dst != end; ++dst, ++src)
				*dst = *src;
		}
        void fill(const ImageData& im)
        {
            Iterator_const& src = im.begin();
            Iterator& dst = this->begin();
            Iterator& end = this->end();
            for (; dst != end; ++dst, ++src)
				*dst = *src;
        }
        /// Write image to file
        virtual void write(const std::string &filename) const = 0;
        /// Get geometrical info
        std::shared_ptr<const VoxelisedGeometricalInfo3D > get_geom_info() const
        {
            // If the geometrical info has not been created yet, throw an error
            if (!_geom_info_sptr) {
                std::cout << "\nGeometrical info not initialised. This implies that your constructor did not call set_up_geom_info().\n";
                throw std::runtime_error("Geometrical info not initialised. This implies that"
                                         " your constructor did not call set_up_geom_info().");
            }
            return _geom_info_sptr;
        }
        /// Clone and return as unique pointer.
        std::unique_ptr<ImageData> clone() const
        {
            return std::unique_ptr<ImageData>(this->clone_impl());
        }
    protected:
        /// Clone helper function. Don't use.
        virtual ImageData* clone_impl() const = 0;
        /// Populate the geometrical info metadata (from the image's own metadata)
        virtual void set_up_geom_info() = 0;
        std::shared_ptr<VoxelisedGeometricalInfo3D> _geom_info_sptr;
	};
}

#endif
