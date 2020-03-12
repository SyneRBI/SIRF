/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

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

/*!
\file
\ingroup STIR Extensions
\brief Specification file for data handling types not present in STIR.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef STIR_DATA_CONTAINER_TYPES
#define STIR_DATA_CONTAINER_TYPES

#include <stdlib.h>

#include <chrono>
#include <fstream>
#include <exception>

#include "sirf/iUtilities/LocalisedException.h"
#include "sirf/iUtilities/DataHandle.h"
#include "sirf/common/DataContainer.h"
#include "sirf/common/ANumRef.h"
#include "sirf/common/PETImageData.h"
#include "sirf/STIR/stir_types.h"
#include "sirf/common/GeometricalInfo.h"
#include "stir/ZoomOptions.h"

namespace sirf {

	class SIRFUtilities {
	public:
		static long long milliseconds()
		{
			auto now = std::chrono::system_clock::now();
			auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
			return (long long)ms.count();
		}
		static std::string scratch_file_name()
		{
			static int calls = 0;
			char buff[32];
			long long int ms = milliseconds();
			calls++;
			sprintf(buff, "tmp_%d_%lld", calls, ms);
			return std::string(buff);
		}
	};

	/*!
	\ingroup STIR Extensions
	\brief STIR ProjDataInterfile wrapper with additional file managing features.

	This derived class has additional capability of deleting the file it handles
	when an object of this class goes out of existence.
	*/

	class ProjDataFile : public stir::ProjDataInterfile {

	public:
		ProjDataFile(const stir::ProjData& pd, const std::string& filename, bool owns_file = true) :
			stir::ProjDataInterfile(pd.get_exam_info_sptr(),
			pd.get_proj_data_info_sptr(),
			filename, std::ios::in | std::ios::out | std::ios::trunc),
			_filename(filename),
			_owns_file(owns_file)
		{}
		ProjDataFile(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info,
			const std::string& filename, bool owns_file = true) :
			stir::ProjDataInterfile(sptr_exam_info, sptr_proj_data_info,
			filename, std::ios::in | std::ios::out | std::ios::trunc),
			_filename(filename),
			_owns_file(owns_file)
		{}
		~ProjDataFile()
		{
			close_stream();
			clear_stream();
			if (!_owns_file)
				return;
			int err;
			err = std::remove((_filename + ".hs").c_str());
			if (err)
				std::cout << "deleting " << _filename << ".hs "
				<< "failed, please delete manually" << std::endl;
			err = std::remove((_filename + ".s").c_str());
			if (err)
				std::cout << "deleting " << _filename << ".s "
				<< "failed, please delete manually" << std::endl;
		}
		stir::shared_ptr<std::iostream> sino_stream_sptr()
		{
			return sino_stream;
		}
		void close_stream()
		{
			((std::fstream*)sino_stream.get())->close();
		}
		void clear_stream()
		{
			((std::fstream*)sino_stream.get())->clear();
		}
	private:
		bool _owns_file;
		std::string _filename;
	};

	/*!
	\ingroup STIR Extensions
	\brief STIR ProjData wrapper with added functionality.

	This class enjoys some features of STIR ProjData and, additionally,
	implements the linear algebra functionality specified by the
	abstract base class aDatacontainer, and provides means for the data
	storage mode (file/memory) selection.
	*/

	class PETAcquisitionData : public DataContainer {
	public:
		virtual ~PETAcquisitionData() {}

		// virtual constructors
		virtual PETAcquisitionData* same_acquisition_data
			(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info) const = 0;
		virtual stir::shared_ptr<PETAcquisitionData> new_acquisition_data() const = 0;

		stir::shared_ptr<PETAcquisitionData> single_slice_rebinned_data(
			const int num_segments_to_combine,
			const int num_views_to_combine = 1,
			const int num_tang_poss_to_trim = 0,
			const bool do_normalisation = true,
			const int max_in_segment_num_to_process = -1
			)
		{
			stir::shared_ptr<stir::ProjDataInfo> out_proj_data_info_sptr(
				stir::SSRB(*data()->get_proj_data_info_ptr(),
				num_segments_to_combine,
				num_views_to_combine,
				num_tang_poss_to_trim,
				max_in_segment_num_to_process
				));
			stir::shared_ptr<PETAcquisitionData> 
				sptr(same_acquisition_data
				(data()->get_exam_info_sptr(), out_proj_data_info_sptr));
			SSRB(*sptr, *data(), do_normalisation);
			return sptr;
		}

		static std::string storage_scheme()
		{
			static bool initialized = false;
			if (!initialized) {
				_storage_scheme = "file";
				initialized = true;
			}
			return _storage_scheme;
		}
		static stir::shared_ptr<PETAcquisitionData> storage_template()
		{
			return _template;
		}

		stir::shared_ptr<stir::ProjData> data()
		{
			return _data;
		}
		const stir::shared_ptr<stir::ProjData> data() const
		{
			return _data;
		}
		void set_data(stir::shared_ptr<stir::ProjData> data)
		{
			_data = data;
		}

		// data import/export
		void fill(float v) { data()->fill(v); }
		void fill(const PETAcquisitionData& ad)
		{
			stir::shared_ptr<stir::ProjData> sptr = ad.data();
			data()->fill(*sptr);
		}
		void fill_from(const float* d) { data()->fill_from(d); }
		void copy_to(float* d) { data()->copy_to(d); }
		std::unique_ptr<PETAcquisitionData> clone() const
		{
			return std::unique_ptr<PETAcquisitionData>(clone_impl());
		}

		// data container methods
		unsigned int items() const { return 1; }
		virtual float norm() const;
		virtual void dot(const DataContainer& a_x, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		virtual void multiply
			(const DataContainer& x, const DataContainer& y);
		virtual void divide
			(const DataContainer& x, const DataContainer& y);
		virtual void inv(float a, const DataContainer& x);
		virtual void write(const std::string &filename) const
		{
			ProjDataFile pd(*data(), filename.c_str(), false);
			pd.fill(*data());
		}

		// ProjData methods
		int get_num_tangential_poss()
		{
			return data()->get_num_tangential_poss();
		}
		int get_num_views()
		{
			return data()->get_num_views();
		}
		int get_num_sinograms()
		{
			return data()->get_num_sinograms();
		}
		int get_num_TOF_bins()
		{
			return 1;
		}
		int get_max_segment_num() const
		{
			return data()->get_max_segment_num();
		}
		stir::SegmentBySinogram<float>
			get_segment_by_sinogram(const int segment_num) const
		{
			return data()->get_segment_by_sinogram(segment_num);
		}
		stir::SegmentBySinogram<float>
			get_empty_segment_by_sinogram(const int segment_num) const
		{
			return data()->get_empty_segment_by_sinogram(segment_num);
		}
		virtual stir::Succeeded set_segment(const stir::SegmentBySinogram<float>& s)
		{
			return data()->set_segment(s);
		}
		stir::shared_ptr<stir::ExamInfo> get_exam_info_sptr() const
		{
			return data()->get_exam_info_sptr();
		}
		stir::shared_ptr<stir::ProjDataInfo> get_proj_data_info_sptr() const
		{
			return data()->get_proj_data_info_sptr();
		}

		// ProjData casts
		operator stir::ProjData&() { return *data(); }
		operator const stir::ProjData&() const { return *data(); }
		operator stir::shared_ptr<stir::ProjData>() { return data(); }

		static stir::shared_ptr<stir::ProjDataInfo>
			proj_data_info_from_scanner(std::string scanner_name,
			int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
		{
			stir::shared_ptr<stir::Scanner> 
				sptr_s(stir::Scanner::get_scanner_from_name(scanner_name));
			//std::cout << "scanner: " << sptr_s->get_name().c_str() << '\n';
			if (boost::iequals(sptr_s->get_name(), "unknown")) {
				throw LocalisedException("Unknown scanner", __FILE__, __LINE__);
			}
			int num_views = sptr_s->get_num_detectors_per_ring() / 2 / view_mash_factor;
			int num_tang_pos = sptr_s->get_max_num_non_arccorrected_bins();
			if (max_ring_diff < 0)
				max_ring_diff = sptr_s->get_num_rings() - 1;
			return std::move(stir::ProjDataInfo::construct_proj_data_info
				(sptr_s, span, max_ring_diff, num_views, num_tang_pos, false));
		}

	protected:
		static std::string _storage_scheme;
		static stir::shared_ptr<PETAcquisitionData> _template;
		stir::shared_ptr<stir::ProjData> _data;
		virtual PETAcquisitionData* clone_impl() const = 0;
		PETAcquisitionData* clone_base() const
		{
			stir::shared_ptr<stir::ExamInfo> sptr_ei = get_exam_info_sptr();
			stir::shared_ptr<stir::ProjDataInfo> sptr_pdi = get_proj_data_info_sptr();
			PETAcquisitionData* ptr = 
				_template->same_acquisition_data(sptr_ei, sptr_pdi);
			ptr->fill(*this);
			return ptr;
		}

	};

	/*!
	\ingroup STIR Extensions
	\brief In-file implementation of PETAcquisitionData.

	*/

	class PETAcquisitionDataInFile : public PETAcquisitionData {
	public:
		PETAcquisitionDataInFile() : _owns_file(false) {}
		PETAcquisitionDataInFile(const char* filename) : _owns_file(false)
		{
			_data = stir::ProjData::read_from_file(filename);
		}
		PETAcquisitionDataInFile(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info)
		{
			_data.reset(new ProjDataFile
				(sptr_exam_info, sptr_proj_data_info,
				_filename = SIRFUtilities::scratch_file_name()));
		}
		PETAcquisitionDataInFile(const stir::ProjData& pd) : _owns_file(true)
		{
			_data.reset(new ProjDataFile
				(pd, _filename = SIRFUtilities::scratch_file_name()));
		}
		PETAcquisitionDataInFile
			(stir::shared_ptr<stir::ExamInfo> sptr_ei, std::string scanner_name,
			int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
		{
			stir::shared_ptr<stir::ProjDataInfo> sptr_pdi =
				PETAcquisitionData::proj_data_info_from_scanner
				(scanner_name, span, max_ring_diff, view_mash_factor);
			ProjDataFile* ptr = new ProjDataFile(sptr_ei, sptr_pdi,
				_filename = SIRFUtilities::scratch_file_name());
			ptr->fill(0.0f);
			_data.reset(ptr);
		}
		stir::shared_ptr<PETAcquisitionData> new_acquisition_data(std::string filename)
		{
			stir::shared_ptr<PETAcquisitionDataInFile> sptr_ad(new PETAcquisitionDataInFile);
			sptr_ad->_data.reset(new ProjDataFile(*data(), filename, false));
			return sptr_ad;
		}

		static void init() {
			static bool initialized = false;
			if (!initialized) {
				_storage_scheme = "file";
				_template.reset(new PETAcquisitionDataInFile());
				initialized = true;
				PETAcquisitionData::storage_scheme();
			}
		}
		static void set_as_template()
		{
			init();
			_storage_scheme = "file";
			_template.reset(new PETAcquisitionDataInFile);
		}

		virtual PETAcquisitionData* same_acquisition_data
			(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info) const
		{
			PETAcquisitionData* ptr_ad =
				new PETAcquisitionDataInFile(sptr_exam_info, sptr_proj_data_info);
			return ptr_ad;
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			init();
			DataContainer* ptr = _template->same_acquisition_data(this->get_exam_info_sptr(),
				this->get_proj_data_info_sptr());
			return new ObjectHandle<DataContainer>
				(stir::shared_ptr<DataContainer>(ptr));
		}
		virtual stir::shared_ptr<PETAcquisitionData> new_acquisition_data() const
		{
			init();
			return stir::shared_ptr < PETAcquisitionData >
				(_template->same_acquisition_data(this->get_exam_info_sptr(),
				this->get_proj_data_info_sptr()));
		}

	private:
		bool _owns_file;
		std::string _filename;
		virtual PETAcquisitionDataInFile* clone_impl() const
		{
			init();
			return (PETAcquisitionDataInFile*)clone_base();
		}
	};

	/*!
	\ingroup STIR Extensions
	\brief In-memory implementation of PETAcquisitionData.

	*/

	class PETAcquisitionDataInMemory : public PETAcquisitionData {
	public:
		PETAcquisitionDataInMemory() {}
		PETAcquisitionDataInMemory(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info)
		{
			_data = stir::shared_ptr<stir::ProjData>
				(new stir::ProjDataInMemory(sptr_exam_info, sptr_proj_data_info));
		}
		PETAcquisitionDataInMemory(const stir::ProjData& pd)
		{
			_data = stir::shared_ptr<stir::ProjData>
				(new stir::ProjDataInMemory(pd.get_exam_info_sptr(),
				pd.get_proj_data_info_sptr()));
		}
		PETAcquisitionDataInMemory
			(stir::shared_ptr<stir::ExamInfo> sptr_ei, std::string scanner_name,
			int span = 1, int max_ring_diff = -1, int view_mash_factor = 1)
		{
			stir::shared_ptr<stir::ProjDataInfo> sptr_pdi =
				PETAcquisitionData::proj_data_info_from_scanner
				(scanner_name, span, max_ring_diff, view_mash_factor);
			stir::ProjDataInMemory* ptr = new stir::ProjDataInMemory(sptr_ei, sptr_pdi);
			ptr->fill(0.0f);
			_data.reset(ptr);
		}

		static void init() 
		{ 
			PETAcquisitionDataInFile::init(); 
		}
		static void set_as_template()
		{
			init();
			_storage_scheme = "memory";
			_template.reset(new PETAcquisitionDataInMemory);
		}

		virtual PETAcquisitionData* same_acquisition_data
			(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info) const
		{
			PETAcquisitionData* ptr_ad =
				new PETAcquisitionDataInMemory(sptr_exam_info, sptr_proj_data_info);
			return ptr_ad;
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			init();
			DataContainer* ptr = _template->same_acquisition_data
				(this->get_exam_info_sptr(), this->get_proj_data_info_sptr());
			return new ObjectHandle<DataContainer>
				(stir::shared_ptr<DataContainer>(ptr));
		}
		virtual stir::shared_ptr<PETAcquisitionData> new_acquisition_data() const
		{
			init();
			return stir::shared_ptr < PETAcquisitionData >
				(_template->same_acquisition_data
				(this->get_exam_info_sptr(), this->get_proj_data_info_sptr()));
		}
	private:
		virtual PETAcquisitionDataInMemory* clone_impl() const
		{
			init();
			return (PETAcquisitionDataInMemory*)clone_base();
		}
	};

	/*!
	\ingroup STIR Extensions
	\brief STIR DiscretisedDensity<3, float> wrapper with added functionality.

	This class enjoys some features of STIR DiscretisedDensity<3, float> and,
	additioanally, implements the linear algebra functionality specified by the
	abstract base class aDatacontainer.
	*/

	typedef Image3DF::full_iterator Image3DFIterator;
	typedef Image3DF::const_full_iterator Image3DFIterator_const;

	//class STIRImageData : public aDataContainer < float > {
	class STIRImageData : public PETImageData { //<Iterator, Iterator_const> {
	public:
		//typedef PETImageData<Iterator, Iterator_const>::Iter BaseIter;
		//typedef PETImageData<Iterator, Iterator_const>::Iter_const BaseIter_const;
		typedef ImageData::Iterator BaseIter;
		typedef ImageData::Iterator_const BaseIter_const;
		class Iterator : public BaseIter {
		public:
			Iterator(const Image3DFIterator& iter) : _iter(iter)
			{}
			Iterator& operator=(const Iterator& iter)
			{
				_iter = iter._iter;
				_ref.copy(iter._ref);
				_sptr_iter = iter._sptr_iter;
				return *this;
			}
			virtual Iterator& operator++()
			{
				++_iter;
				return *this;
			}
			//virtual Iterator& operator++(int)
			//{
			//	_sptr_iter.reset(new Iterator(_iter));
			//	++_iter;
			//	return *_sptr_iter;
			//}
			virtual bool operator==(const BaseIter& an_iter) const
			{
				const Iterator& iter = (const Iterator&)an_iter;
				return _iter == iter._iter;
			}
			virtual bool operator!=(const BaseIter& an_iter) const
			{
				const Iterator& iter = (const Iterator&)an_iter;
				return _iter != iter._iter;
			}
			virtual FloatRef& operator*()
			{
				float& v = *_iter;
				_ref.set_ptr(&v);
				return _ref;
			}
		private:
			Image3DFIterator _iter;
			FloatRef _ref;
			std::shared_ptr<Iterator> _sptr_iter;
		};
		class Iterator_const : public BaseIter_const {
		public:
			Iterator_const(const Image3DFIterator_const& iter) : _iter(iter)
			{}
			Iterator_const& operator=(const Iterator_const& iter)
			{
				_iter = iter._iter;
				_ref.copy(iter._ref);
				_sptr_iter = iter._sptr_iter;
				return *this;
			}
			virtual Iterator_const& operator++()
			{
				++_iter;
				return *this;
			}
			//virtual Iterator_const& operator++(int)
			//{
			//	_sptr_iter.reset(new Iterator_const(_iter));
			//	++_iter;
			//	return *_sptr_iter;
			//}
			virtual bool operator==(const BaseIter_const& an_iter) const
			{
				const Iterator_const& iter = (const Iterator_const&)an_iter;
				return _iter == iter._iter;
			}
			virtual bool operator!=(const BaseIter_const& an_iter) const
			{
				const Iterator_const& iter = (const Iterator_const&)an_iter;
				return _iter != iter._iter;
			}
			virtual const FloatRef& operator*() const
			{
				const float& v = *_iter;
				_ref.set_ptr((void*)&v);
				return _ref;
			}
		private:
			Image3DFIterator_const _iter;
			mutable FloatRef _ref;
			std::shared_ptr<Iterator_const> _sptr_iter;
		};
		STIRImageData() {}
        STIRImageData(const ImageData& id);
		STIRImageData(const STIRImageData& image)
		{
			_data.reset(image.data().clone());
            this->set_up_geom_info();
		}
		STIRImageData(const PETAcquisitionData& ad)
		{
			_data.reset(new Voxels3DF(ad.get_exam_info_sptr(),*ad.get_proj_data_info_sptr()));
            this->set_up_geom_info();
		}
		STIRImageData(const Image3DF& image)
		{
			_data.reset(image.clone());
            this->set_up_geom_info();
		}
		STIRImageData(const Voxels3DF& v)
		{
			_data.reset(v.clone());
            this->set_up_geom_info();
		}
		STIRImageData(const stir::ProjDataInfo& pdi)
		{
			_data.reset(new Voxels3DF(pdi));
            this->set_up_geom_info();
		}
		STIRImageData(stir::shared_ptr<Image3DF> ptr)
		{
			_data = ptr;
            this->set_up_geom_info();
		}
		STIRImageData(std::string filename)
		{
			_data = stir::read_from_file<Image3DF>(filename);
            this->set_up_geom_info();
		}
		STIRImageData* same_image_data() const
		{
			STIRImageData* ptr_image = new STIRImageData;
			ptr_image->_data.reset(_data->get_empty_copy());
            ptr_image->set_up_geom_info();
			return ptr_image;
		}
		stir::shared_ptr<STIRImageData> new_image_data()
		{
			return stir::shared_ptr<STIRImageData>(same_image_data());
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(stir::shared_ptr<DataContainer>(same_image_data()));
		}
		unsigned int items() const
		{
			return 1;
		}
        /// Write to file
        virtual void write(const std::string &filename) const;
        /// Write to file using format file.
        /*! This allows speciyfing the output file format used by STIR using a text file.
            Keywords are specified by STIR.

            If "" is passed as argument for format_file, the default format will be used.

            An example is given below for writing the image in the nifti format. STIR uses
            ITK to do this, so ensure that STIR is built with ITK if you wish to use it.
            \verbatim
            OutputFileFormat Parameters:=
                output file format type := ITK
                ITK Output File Format Parameters:=
                number format := float
                number_of_bytes_per_pixel:=4
                default extension:=.nii
                End ITK Output File Format Parameters:=
           End:=
           \endverbatim
        */
        virtual void write(const std::string &filename, const std::string &format_file) const;

		virtual float norm() const;
		virtual void dot(const DataContainer& a_x, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		virtual void multiply(const DataContainer& x,
			const DataContainer& y);
		virtual void divide(const DataContainer& x,
			const DataContainer& y);

		Image3DF& data()
		{
			return *_data;
		}
		const Image3DF& data() const
		{
			return *_data;
		}
		Image3DF* data_ptr()
		{
			return _data.get();
		}
		const Image3DF* data_ptr() const
		{
			return _data.get();
		}
		stir::shared_ptr<Image3DF> data_sptr()
		{
			return _data;
		}
		void set_data_sptr(stir::shared_ptr<Image3DF> sptr_data)
		{
			_data = sptr_data;
		}
		void fill(float v)
		{
			_data->fill(v);
		}
		virtual Dimensions dimensions() const
		{
			Dimensions dim;
			int d[4];
			get_dimensions(d);
			dim["z"] = d[0];
			dim["y"] = d[1];
			dim["x"] = d[2];
			return dim;
		}
		int get_dimensions(int* dim) const;
		void get_voxel_sizes(float* vsizes) const;
		virtual void get_data(float* data) const;
		virtual void set_data(const float* data);
		virtual Iterator& begin()
		{
			_begin.reset(new Iterator(data().begin_all()));
			return *_begin;
		}
		virtual Iterator_const& begin() const
		{
			_begin_const.reset(new Iterator_const(data().begin_all()));
			return *_begin_const;
		}
		virtual Iterator& end()
		{
			_end.reset(new Iterator(data().end_all()));
			return *_end;
		}
		virtual Iterator_const& end() const
		{
			_end_const.reset(new Iterator_const(data().end_all()));
			return *_end_const;
		}

        /// Clone and return as unique pointer.
        std::unique_ptr<STIRImageData> clone() const
        {
            return std::unique_ptr<STIRImageData>(this->clone_impl());
        }

        /// Zoom the image (modifies itself).
        /// All indices and coordinates should be (z,y,x) order.
        /// To leave the size unchanged in any dimension, set the corresponding size to -1.
        void zoom_image(const Coord3DF &zooms={1.f,1.f,1.f}, const Coord3DF &offsets_in_mm={0.f,0.f,0.f},
                        const Coord3DI &new_sizes={-1,-1,-1}, const char * const zoom_options_str="preserve_sum");

        /// Zoom the image (modifies itself).
        /// All indices and coordinates should be (z,y,x) order.
        /// To leave the size unchanged in any dimension, set the corresponding size to -1.
        void zoom_image(const Coord3DF &zooms={1.f,1.f,1.f}, const Coord3DF &offsets_in_mm={0.f,0.f,0.f},
                        const Coord3DI &new_sizes={-1,-1,-1},
                        const stir::ZoomOptions zoom_options=stir::ZoomOptions::preserve_sum);

        /// Move to scanner centre. The acquisition needs to be supplied such that in the future,
        /// bed offset etc can be taken into account.
        void move_to_scanner_centre(const PETAcquisitionData &);

        /// Populate the geometrical info metadata (from the image's own metadata)
        virtual void set_up_geom_info();

    private:
        /// Clone helper function. Don't use.
        virtual STIRImageData* clone_impl() const
        {
            return new STIRImageData(*this);
        }

	protected:

		stir::shared_ptr<Image3DF> _data;
		mutable stir::shared_ptr<Iterator> _begin;
		mutable stir::shared_ptr<Iterator> _end;
		mutable stir::shared_ptr<Iterator_const> _begin_const;
		mutable stir::shared_ptr<Iterator_const> _end_const;
	};

}  // namespace sirf

#endif
