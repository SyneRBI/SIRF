/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 Rutherford Appleton Laboratory STFC
Copyright 2018 University College London

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/
#pragma once

#ifndef SIRF_ABSTRACT_IMAGE_DATA_TYPE
#define SIRF_ABSTRACT_IMAGE_DATA_TYPE

#include "sirf/common/ANumRef.h"
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
        //virtual void write(const std::string &filename) const = 0;
		bool operator==(const ImageData& id) const
		{
			GeometricalInfo<3, 3>& gi_self = (GeometricalInfo<3, 3>&)*get_geom_info_sptr();
			GeometricalInfo<3, 3>& gi_other = (GeometricalInfo<3, 3>&)*id.get_geom_info_sptr();
			if (gi_self != gi_other)
				return false;
			float s = 0.0f;
			float sx = 0.0f;
			float sy = 0.0f;
			complex_float_t zx;
			complex_float_t zy;
			Iterator_const& x = this->begin();
			Iterator_const& y = id.begin();
			for (; x != this->end(); ++x, ++y) {
				zx = (*x).complex_float();
				zy = (*y).complex_float();
				sx += std::abs(zx*zx);
				sy += std::abs(zy*zy);
				zx -= zy;
				s += std::abs(zx*zx);
			}
			float t = std::max(sx, sy);
			bool same = (s <= 1e-6*t);
			return same;
		}
		bool operator!=(const ImageData& id) const
		{
			return !(*this == id);
		}
		/// Get geometrical info
        std::shared_ptr<const VoxelisedGeometricalInfo3D > get_geom_info_sptr() const
        {
            // If the geometrical info has not been created yet, throw an error
            if (!_geom_info_sptr)
                throw std::runtime_error("Geometrical info not initialised. This implies that"
                                         " your constructor did not call set_up_geom_info() or there was "
                                         "an error. Build in debug mode and lookout for any printed text "
                                         "containing '::set_up_geom_info()'.");
            return _geom_info_sptr;
        }
        /// Clone and return as unique pointer.
        std::unique_ptr<ImageData> clone() const
        {
            return std::unique_ptr<ImageData>(this->clone_impl());
        }
        /// Is complex? Unless overwridden (Gadgetron), assume not complex.
        virtual bool is_complex() const { return false; }
        /// Reorient image. Requires that dimesions and spacing match
        virtual void reorient(const VoxelisedGeometricalInfo3D &);
        /// Can reorient? (check dimensions and spacing)
        static bool can_reorient(const VoxelisedGeometricalInfo3D &geom_1, const VoxelisedGeometricalInfo3D &geom_2, const bool throw_error);
        /// Populate the geometrical info metadata (from the image's own metadata)
        virtual void set_up_geom_info() = 0;
    protected:
        /// Clone helper function. Don't use.
        virtual ImageData* clone_impl() const = 0;
        /// Set geom info
        void set_geom_info(const std::shared_ptr<VoxelisedGeometricalInfo3D> geom_info_sptr) { _geom_info_sptr = geom_info_sptr; }
    private:
        std::shared_ptr<VoxelisedGeometricalInfo3D> _geom_info_sptr;
	};
}

#endif