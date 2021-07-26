/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
Copyright 2018 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

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
\ingroup PET
\brief Specification file for data handling types not present in STIR.

\author Evgueni Ovtchinnikov
\author Richard Brown
\author SyneRBI
*/

#ifndef STIR_DATA_CONTAINER_TYPES
#define STIR_DATA_CONTAINER_TYPES

#include <stdlib.h>

#include <chrono>
#include <fstream>
#include <exception>

#include "sirf/iUtilities/LocalisedException.h"
#include "sirf/iUtilities/DataHandle.h"
#include "sirf/common/iequals.h"
#include "sirf/common/JacobiCG.h"
#include "sirf/common/DataContainer.h"
#include "sirf/common/ANumRef.h"
#include "sirf/common/PETImageData.h"
#include "sirf/STIR/stir_types.h"
#include "sirf/common/GeometricalInfo.h"
#include "stir/ZoomOptions.h"

#if STIR_VERSION < 050000
#define SPTR_WRAP(X) X->create_shared_clone()
#else
#define SPTR_WRAP(X) X
#endif

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
	\ingroup PET
	\brief STIR ProjDataInterfile wrapper with additional file managing features.

	This derived class has additional capability of deleting the file it handles
	when an object of this class goes out of existence.
	*/

	class ProjDataFile : public stir::ProjDataInterfile {

	public:
		ProjDataFile(const stir::ProjData& pd, const std::string& filename, bool owns_file = true) :
			stir::ProjDataInterfile(pd.get_exam_info_sptr(),
			pd.get_proj_data_info_sptr()->create_shared_clone(),
			filename, std::ios::in | std::ios::out | std::ios::trunc),
			_filename(filename),
			_owns_file(owns_file)
		{}
		ProjDataFile(stir::shared_ptr<stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info,
			const std::string& filename, bool owns_file = true) :
			stir::ProjDataInterfile(SPTR_WRAP(sptr_exam_info), SPTR_WRAP(sptr_proj_data_info),
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
	\ingroup PET
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
			(stir::shared_ptr<const stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info) const = 0;
		virtual stir::shared_ptr<PETAcquisitionData> new_acquisition_data() const = 0;

		virtual bool is_complex() const
		{
			return false;
		}

		//! rebin the data to lower resolution by adding
		/*!
		  \param num_segments_to_combine combines multiple oblique 'segments' together. If set to the
		    total number of segments, this corresponds to SSRB. Another example is if the input data
			has 'span=1', the output span will be equal to the \c num_segments_to_combine.
		  \param num_views_to_combine combines neighbouring views. Needs to be a divisor of the total
		    number of views in the data.
		  \param num_tang_poss_to_trim removes a number of tangential positions (horizontal direction
		    in the sinogram) at each end
		  \param do_normalisation if \c true, averages the data, otherwise it adds the data. Often
		    the latter is required for emission data (as it preserves Poisson statistics),
			while the former should be used for corrected data (or for attenuation correction factors).
		  \param max_in_segment_num_to_process by default all input data are used. If set to a non-negative
		    number, it will remove the most oblique segments.
		*/
		stir::shared_ptr<PETAcquisitionData> single_slice_rebinned_data(
			const int num_segments_to_combine,
			const int num_views_to_combine = 1,
			const int num_tang_poss_to_trim = 0,
			const bool do_normalisation = true,
			const int max_in_segment_num_to_process = -1
			)
		{
			stir::shared_ptr<stir::ProjDataInfo> out_proj_data_info_sptr(
				stir::SSRB(*data()->get_proj_data_info_sptr(),
				num_segments_to_combine,
				num_views_to_combine,
				num_tang_poss_to_trim,
				max_in_segment_num_to_process
				));
			stir::shared_ptr<PETAcquisitionData> 
				sptr(same_acquisition_data
                                     (this->get_exam_info_sptr(), out_proj_data_info_sptr));
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
		virtual void fill(const float v) { data()->fill(v); }
		virtual void fill(const PETAcquisitionData& ad)
		{
			if (ad.is_empty())
				THROW("The source of PETAcquisitionData::fill is empty");
			stir::shared_ptr<stir::ProjData> sptr = ad.data();
			data()->fill(*sptr);
		}
		virtual void fill_from(const float* d) { data()->fill_from(d); }
		virtual void copy_to(float* d) const { data()->copy_to(d); }
		std::unique_ptr<PETAcquisitionData> clone() const
		{
			return std::unique_ptr<PETAcquisitionData>(clone_impl());
		}

		// data container methods
		unsigned int items() const {
			if (_is_empty != -1)
				return _is_empty ? 0 : 1;
			try {
				get_segment_by_sinogram(0);
			}
			catch (std::string msg) {
				_is_empty = 1;
				return 0; // no data found - this must be a template
			}
			_is_empty = 0;
			return 1; // data ok
		}
		virtual float norm() const;
		virtual void dot(const DataContainer& a_x, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		virtual void xapyb(
			const DataContainer& a_x, const void* ptr_a,
			const DataContainer& a_y, const void* ptr_b);	
		virtual void xapyb(
			const DataContainer& a_x, const DataContainer& a_a,
			const DataContainer& a_y, const DataContainer& a_b);
		virtual void multiply(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 1);
		}
		virtual void divide(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 2);
		}
		virtual void maximum(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 3);
		}
		virtual void minimum(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 4);
		}
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
		//! total number of (2D) sinograms
		/*! note that for TOF data, this includes the TOF bins.
		    \see get_num_non_TOF_sinograms()
	    */
		int get_num_sinograms() const
		{
			return data()->get_num_sinograms();
        }
		//! total number of (2D) sinograms ignoring time-of-flight
		/*! This does include the oblique data as well. */
		int get_num_non_TOF_sinograms() const
		{
			return data()->get_num_non_tof_sinograms();
		}

		int get_num_TOF_bins()
		{
			return data()->get_num_tof_poss();
		}
		size_t get_dimensions(int* dim)
		{
			dim[0] = get_num_tangential_poss();
			dim[1] = get_num_views();
			dim[2] = get_num_non_TOF_sinograms();
			dim[3] = get_num_TOF_bins();
			return static_cast<size_t>(dim[0] * dim[1] * dim[2] * dim[3]);
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
		stir::shared_ptr<const stir::ExamInfo> get_exam_info_sptr() const
		{
			return data()->get_exam_info_sptr();
		}
		stir::shared_ptr<const stir::ProjDataInfo> get_proj_data_info_sptr() const
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
			if (sirf::iequals(sptr_s->get_name(), "unknown")) {
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
			stir::shared_ptr<stir::ProjDataInfo> sptr_pdi = this->get_proj_data_info_sptr()->create_shared_clone();
			PETAcquisitionData* ptr = 
				_template->same_acquisition_data(this->get_exam_info_sptr(), sptr_pdi);
			if (!this->is_empty())
				ptr->fill(*this);
			return ptr;
		}

	private:
		mutable int _is_empty = -1;
		void binary_op_(const DataContainer& a_x, const DataContainer& a_y, int job);
	};

	/*!
	\ingroup PET
	\brief In-file implementation of PETAcquisitionData.

	*/

	class PETAcquisitionDataInFile : public PETAcquisitionData {
	public:
		PETAcquisitionDataInFile() : _owns_file(false) {}
		PETAcquisitionDataInFile(const char* filename) : _owns_file(false)
		{
			_data = stir::ProjData::read_from_file(filename);
		}
		PETAcquisitionDataInFile(stir::shared_ptr<const stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info)
		{
			_data.reset(new ProjDataFile
                                    (MAKE_SHARED<stir::ExamInfo>(*sptr_exam_info), sptr_proj_data_info,
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
			(stir::shared_ptr<const stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<stir::ProjDataInfo> sptr_proj_data_info) const
		{
			PETAcquisitionData* ptr_ad =
				new PETAcquisitionDataInFile(sptr_exam_info, sptr_proj_data_info);
			return ptr_ad;
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			init();
			DataContainer* ptr = _template->same_acquisition_data(
                                this->get_exam_info_sptr(),
				this->get_proj_data_info_sptr()->create_shared_clone());
			return new ObjectHandle<DataContainer>
				(stir::shared_ptr<DataContainer>(ptr));
		}
		virtual stir::shared_ptr<PETAcquisitionData> new_acquisition_data() const
		{
			init();
			return stir::shared_ptr < PETAcquisitionData >
				(_template->same_acquisition_data(this->get_exam_info_sptr(),
				this->get_proj_data_info_sptr()->create_shared_clone()));
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
	\ingroup PET
	\brief In-memory implementation of PETAcquisitionData.

	*/

	class PETAcquisitionDataInMemory : public PETAcquisitionData {
	public:
		PETAcquisitionDataInMemory() {}
		PETAcquisitionDataInMemory(stir::shared_ptr<const stir::ExamInfo> sptr_exam_info,
			stir::shared_ptr<const stir::ProjDataInfo> sptr_proj_data_info)
		{
			_data = stir::shared_ptr<stir::ProjData>
				(new stir::ProjDataInMemory(SPTR_WRAP(sptr_exam_info), SPTR_WRAP(sptr_proj_data_info)));
		}
		PETAcquisitionDataInMemory(const stir::ProjData& templ)
		{
			_data = stir::shared_ptr<stir::ProjData>
				(new stir::ProjDataInMemory(templ.get_exam_info_sptr(),
					templ.get_proj_data_info_sptr()->create_shared_clone()));
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
        /// Constructor for PETAcquisitionDataInMemory from filename
        PETAcquisitionDataInMemory(const char* filename)
        {
            auto pd_sptr = stir::ProjData::read_from_file(filename);
			bool is_empty = false;
			try {
				pd_sptr->get_segment_by_sinogram(0);
			}
			catch (...) {
				is_empty = true;
			}
			if (is_empty)
				_data = stir::shared_ptr<stir::ProjData>
					(new stir::ProjDataInMemory(pd_sptr->get_exam_info_sptr(),
						pd_sptr->get_proj_data_info_sptr()->create_shared_clone()));
			else
				_data = stir::shared_ptr<stir::ProjData>
				(new stir::ProjDataInMemory(*pd_sptr));
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
			(stir::shared_ptr<const stir::ExamInfo> sptr_exam_info,
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
				(this->get_exam_info_sptr(),
                                 this->get_proj_data_info_sptr()->create_shared_clone());
			return new ObjectHandle<DataContainer>
				(stir::shared_ptr<DataContainer>(ptr));
		}
		virtual stir::shared_ptr<PETAcquisitionData> new_acquisition_data() const
		{
			init();
			return stir::shared_ptr < PETAcquisitionData >
				(_template->same_acquisition_data
				(this->get_exam_info_sptr(),
                                 this->get_proj_data_info_sptr()->create_shared_clone()));
		}
        /// fill with single value
        virtual void fill(const float v)
        {
            stir::ProjDataInMemory *pd_ptr = dynamic_cast<stir::ProjDataInMemory*>(data().get());
            // If cast failed, fall back to general method
            if (is_null_ptr(pd_ptr))
                return this->PETAcquisitionData::fill(v);

            // do it
            auto iter = pd_ptr->begin();
            while (iter != pd_ptr->end())
                *iter++ = v;
        }
        /// fill from another PETAcquisitionData
        virtual void fill(const PETAcquisitionData& ad)
        {
            // Can only do this if both are PETAcquisitionDataInMemory
            stir::ProjDataInMemory *pd_ptr = dynamic_cast<stir::ProjDataInMemory*>(data().get());
            const stir::ProjDataInMemory *pd2_ptr = dynamic_cast<const stir::ProjDataInMemory*>(ad.data().get());
            // If either cast failed, fall back to general method
            if (is_null_ptr(pd_ptr) || is_null_ptr(pd2_ptr))
                return this->PETAcquisitionData::fill(ad);

            // do it
            auto iter = pd_ptr->begin();
            auto iter_other = pd2_ptr->begin();
            while (iter != pd_ptr->end())
                *iter++ = *iter_other++;
        }
        /// Fill from float array
        virtual void fill_from(const float* d)
        {
            stir::ProjDataInMemory *pd_ptr = dynamic_cast<stir::ProjDataInMemory*>(data().get());
            // If cast failed, fall back to general method
            if (is_null_ptr(pd_ptr))
                return this->PETAcquisitionData::fill_from(d);

            // do it
            auto iter = pd_ptr->begin();
            while (iter != pd_ptr->end())
                *iter++ = *d++;
        }
        /// Copy to float array
        virtual void copy_to(float* d) const
        {
            const stir::ProjDataInMemory *pd_ptr = dynamic_cast<const stir::ProjDataInMemory*>(data().get());
            // If cast failed, fall back to general method
            if (is_null_ptr(pd_ptr))
                return this->PETAcquisitionData::copy_to(d);

            // do it
            auto iter = pd_ptr->begin();
            while (iter != pd_ptr->end())
                *d++ = *iter++;
        }
        /// Calculate the norm
        virtual float norm() const
        {
            const stir::ProjDataInMemory *pd_ptr = dynamic_cast<const stir::ProjDataInMemory*>(data().get());
            // If cast failed, fall back to general method
            if (is_null_ptr(pd_ptr))
                return this->PETAcquisitionData::norm();

            // do it
            double t = 0.0;
            auto iter = pd_ptr->begin();
            while (iter != pd_ptr->end())
                t += (*iter) * (*iter++);
            return sqrt((float)t);
        }
        /// Dot between "this" and "other"
        virtual void dot(const DataContainer& a_x, void* ptr) const
        {
            auto x = dynamic_cast<const PETAcquisitionData*>(&a_x);
            // Can only do this if both are PETAcquisitionDataInMemory
            stir::ProjDataInMemory *pd_ptr = dynamic_cast<stir::ProjDataInMemory*>(data().get());
            const stir::ProjDataInMemory *pd2_ptr = dynamic_cast<const stir::ProjDataInMemory*>(x->data().get());
            // If either cast failed, fall back to general method
            if (is_null_ptr(pd_ptr) || is_null_ptr(pd2_ptr))
                return this->PETAcquisitionData::dot(a_x,ptr);

            // do it
            double t = 0.0;
            auto iter = pd_ptr->begin();
            auto iter_other = pd2_ptr->begin();
            while (iter != pd_ptr->end())
                t += (*iter++) * (*iter_other++);

            float* ptr_t = (float*)ptr;
            *ptr_t = (float)t;
        }
        /// Element-wise multiplication of x and y. Store result in "this"
        virtual void multiply(const DataContainer& x, const DataContainer& y)
        {
            auto a_x = dynamic_cast<const PETAcquisitionData*>(&x);
            auto a_y = dynamic_cast<const PETAcquisitionData*>(&y);

            // Can only do this if all are PETAcquisitionDataInMemory
            auto *pd_ptr   = dynamic_cast<stir::ProjDataInMemory*>(data().get());
            auto *pd_x_ptr = dynamic_cast<const stir::ProjDataInMemory*>(a_x->data().get());
            auto *pd_y_ptr = dynamic_cast<const stir::ProjDataInMemory*>(a_y->data().get());

            // If either cast failed, fall back to general method
            if (is_null_ptr(pd_ptr) || is_null_ptr(pd_x_ptr) || is_null_ptr(pd_x_ptr))
                return this->PETAcquisitionData::multiply(x,y);

            // do it
            auto iter = pd_ptr->begin();
            auto iter_x = pd_x_ptr->begin();
            auto iter_y = pd_y_ptr->begin();
            while (iter != pd_ptr->end())
                *iter++ = (*iter_x++) * (*iter_y++);
        }
        /// Element-wise division of x and y. Store result in "this"
        virtual void divide(const DataContainer& x, const DataContainer& y)
        {
            auto a_x = dynamic_cast<const PETAcquisitionData*>(&x);
            auto a_y = dynamic_cast<const PETAcquisitionData*>(&y);

            // Can only do this if all are PETAcquisitionDataInMemory
            auto *pd_ptr   = dynamic_cast<stir::ProjDataInMemory*>(data().get());
            auto *pd_x_ptr = dynamic_cast<const stir::ProjDataInMemory*>(a_x->data().get());
            auto *pd_y_ptr = dynamic_cast<const stir::ProjDataInMemory*>(a_y->data().get());

            // If either cast failed, fall back to general method
            if (is_null_ptr(pd_ptr) || is_null_ptr(pd_x_ptr) || is_null_ptr(pd_x_ptr))
                return this->PETAcquisitionData::divide(x,y);

            // do it
            auto iter = pd_ptr->begin();
            auto iter_x = pd_x_ptr->begin();
            auto iter_y = pd_y_ptr->begin();
            while (iter != pd_ptr->end())
                *iter++ = (*iter_x++) / (*iter_y++);
        }

	private:
		virtual PETAcquisitionDataInMemory* clone_impl() const
		{
			init();
			return (PETAcquisitionDataInMemory*)clone_base();
		}
	};

	/*!
	\ingroup PET
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
                  _data.reset(new Voxels3DF(MAKE_SHARED<stir::ExamInfo>(*ad.get_exam_info_sptr()),*ad.get_proj_data_info_sptr()));
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
		virtual bool is_complex() const
		{
			return false;
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
		virtual void xapyb(
			const DataContainer& a_x, const void* ptr_a,
			const DataContainer& a_y, const void* ptr_b);		
		virtual void xapyb(
			const DataContainer& a_x, const DataContainer& a_a,
			const DataContainer& a_y, const DataContainer& a_b);
		virtual void multiply(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 1);
		}
		virtual void divide(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 2);
		}
		virtual void maximum(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 3);
		}
		virtual void minimum(const DataContainer& x, const DataContainer& y)
		{
			binary_op_(x, y, 4);
		}

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
		stir::shared_ptr<const Image3DF> data_sptr() const
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
		void scale(float s);
		float dot(const DataContainer& a_x) const
		{
			float s;
			dot(a_x, &s);
			return s;
		}
		void axpby(
			float a, const DataContainer& a_x,
			float b, const DataContainer& a_y)
		{	
			axpby(&a, a_x, &b, a_y);
		}
		void xapyb(
			const DataContainer& a_x, float a,
			const DataContainer& a_y, float b)
		{
			xapyb(a_x, &a, a_y, &b);
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
		void binary_op_(const DataContainer& a_x, const DataContainer& a_y, int job);

	protected:

		stir::shared_ptr<Image3DF> _data;
		mutable stir::shared_ptr<Iterator> _begin;
		mutable stir::shared_ptr<Iterator> _end;
		mutable stir::shared_ptr<Iterator_const> _begin_const;
		mutable stir::shared_ptr<Iterator_const> _end_const;
	};

}  // namespace sirf

#endif
