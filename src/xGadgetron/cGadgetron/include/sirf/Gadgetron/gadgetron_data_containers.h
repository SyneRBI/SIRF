/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London
Copyright 2020 Physikalisch-Technische Bundesanstalt (PTB)

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
\ingroup Gadgetron Data Containers
\brief Specification file for data container classes for Gadgetron data.

\author Evgueni Ovtchinnikov
\author Johannes Mayer
\author SyneRBI
*/

#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include "sirf/common/DataContainer.h"
#include "sirf/common/MRImageData.h"
#include "sirf/common/multisort.h"
#include "sirf/Gadgetron/ismrmrd_fftw.h"
#include "sirf/Gadgetron/cgadgetron_shared_ptr.h"
#include "sirf/Gadgetron/gadgetron_image_wrap.h"
#include "sirf/iUtilities/LocalisedException.h"

//#define DYNAMIC_CAST(T, X, Y) T& X = (T&)Y
#define DYNAMIC_CAST(T, X, Y) T& X = dynamic_cast<T&>(Y)

/*!
\ingroup Gadgetron Data Containers
\brief Acquisitions filter.

Some acquisitions do not participate directly in the reconstruction process
(e.g. noise calibration acquisitions).
*/
#define TO_BE_IGNORED(acq) \
	(!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION) && \
	!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING) && \
	!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_MEASUREMENT) && \
	!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_REVERSE) && \
	(acq).flags() >= (1 << (ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT - 1)))

/*!
\ingroup Gadgetron Data Containers
\brief Serialized ISMRMRD acquisition header (cf. ismrmrd.h).

*/

namespace sirf {

	class AcquisitionsInfo {
	public:
		AcquisitionsInfo(std::string data = "") : data_(data)
        {
			if (data.empty())
				have_header_ = false;
			else {
				deserialize();
				have_header_ = true;
			}
        }
		AcquisitionsInfo& operator=(std::string data)
		{
			data_ = data;
			if (data.empty())
				have_header_ = false;
			else {
				deserialize();
				have_header_ = true;
			}
			return *this;
		}
		const char* c_str() const { return data_.c_str(); }
		operator std::string&() { return data_; }
		operator const std::string&() const { return data_; }
        bool empty() const { return data_.empty(); }
        const ISMRMRD::IsmrmrdHeader& get_IsmrmrdHeader() const 
		{
			if (!have_header_)
				deserialize();
			return header_; 
		}

	private:
		void deserialize() const
		{
			if (!this->empty())
				ISMRMRD::deserialize(data_.c_str(), header_);
		}
		std::string data_;
        mutable ISMRMRD::IsmrmrdHeader header_;
		bool have_header_;
	};

    class KSpaceSorting
    {
        static int const num_kspace_dims_ = 7 + ISMRMRD::ISMRMRD_Constants::ISMRMRD_USER_INTS;

    public:

        typedef std::array<int, num_kspace_dims_> TagType;
        typedef std::vector<int> SetType;

        KSpaceSorting(){
            for(int i=0; i<num_kspace_dims_; ++i)
                this->tag_[i] = -1;
        }

        KSpaceSorting(TagType tag){
            this->tag_ = tag;
            this->idx_set_ = {};
        }

        KSpaceSorting(TagType tag, SetType idx_set){
            this->tag_ = tag;
            this->idx_set_ = idx_set;
        }

        TagType get_tag(void) const {return tag_;}
        SetType get_idx_set(void) const {return idx_set_;}
        void add_idx_to_set(size_t const idx){this->idx_set_.push_back(idx);}

        static TagType get_tag_from_acquisition(ISMRMRD::Acquisition acq)
        {
            TagType tag;
            tag[0] = acq.idx().average;
            tag[1] = acq.idx().slice;
            tag[2] = acq.idx().contrast;
            tag[3] = acq.idx().phase;
            tag[4] = acq.idx().repetition;
            tag[5] = acq.idx().set;
            tag[6] = 0; //acq.idx().segment;

            for(int i=7; i<tag.size(); ++i)
                tag[i]=acq.idx().user[i];

            return tag;
        }

        bool is_first_set() const {
            bool is_first= (tag_[0] == 0);
            if(is_first)
            {
               for(int dim=2; dim<num_kspace_dims_; ++dim)
                   is_first *= (tag_[dim] == 0);
            }
            return is_first;
        }

    private:

        // order is [average, slice, contrast, phase, repetition, set, segment, user_ (0,...,ISMRMRD_USER_INTS-1)]
        TagType tag_;
        SetType idx_set_;

    };

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract MR acquisition data container class.

	*/
	class MRAcquisitionData : public DataContainer {
	public:
		// static methods

		static std::string storage_scheme()
		{
			static bool initialized = false;
			if (!initialized) {
				_storage_scheme = "file";
				initialized = true;
			}
			return _storage_scheme;
		}

		// ISMRMRD acquisitions algebra: acquisitions viewed as vectors of 
		// acquisition data
		// y := a x + b y
		static void axpby
			(complex_float_t a, const ISMRMRD::Acquisition& acq_x,
			complex_float_t b, ISMRMRD::Acquisition& acq_y);
		// the inner (l2) product of x and y
		static complex_float_t dot
			(const ISMRMRD::Acquisition& acq_x, const ISMRMRD::Acquisition& acq_y);
		// elementwise multiplication
		// y := x .* y
		static void multiply
			(const ISMRMRD::Acquisition& acq_x, ISMRMRD::Acquisition& acq_y);
		// elementwise division
		// y := x ./ y
		static void divide
			(const ISMRMRD::Acquisition& acq_x, ISMRMRD::Acquisition& acq_y);
		// l2 norm of x
		static float norm(const ISMRMRD::Acquisition& acq_x);

		// abstract methods

		virtual void empty() = 0;
		virtual void take_over(MRAcquisitionData&) = 0;

		// the number of acquisitions in the container
		virtual unsigned int number() const = 0;

		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const = 0;
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
		virtual void append_acquisition(ISMRMRD::Acquisition& acq) = 0;

		virtual void copy_acquisitions_info(const MRAcquisitionData& ac) = 0;
		virtual void copy_acquisitions_data(const MRAcquisitionData& ac) = 0;

		// 'export' constructors: workaround for creating 'ABC' objects
		virtual gadgetron::unique_ptr<MRAcquisitionData> new_acquisitions_container() = 0;
		virtual MRAcquisitionData*
			same_acquisitions_container(const AcquisitionsInfo& info) const = 0;

		virtual void set_data(const complex_float_t* z, int all = 1) = 0;
		virtual void get_data(complex_float_t* z, int all = 1);

		// acquisition data algebra
		virtual void dot(const DataContainer& dc, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		virtual void multiply(
			const DataContainer& a_x,
			const DataContainer& a_y);
		virtual void divide(
			const DataContainer& a_x,
			const DataContainer& a_y);
		virtual float norm() const;

		virtual void write(const std::string &filename) const;

		// regular methods

		AcquisitionsInfo acquisitions_info() const { return acqs_info_; }
		void set_acquisitions_info(std::string info) { acqs_info_ = info; }

		gadgetron::unique_ptr<MRAcquisitionData> clone() const
		{
			return gadgetron::unique_ptr<MRAcquisitionData>(this->clone_impl());
		}

		bool undersampled() const;
		int get_acquisitions_dimensions(size_t ptr_dim) const;
	
		void sort();
		void sort_by_time();
		bool sorted() const { return sorted_; }
		void set_sorted(bool sorted) { sorted_ = sorted; }

        std::vector<std::vector<int> > get_kspace_order(const bool get_first_subset_order=false) const;
        void organise_kspace();

        virtual void get_subset(MRAcquisitionData& subset, const std::vector<int> subset_idx) const;
        virtual void set_subset(const MRAcquisitionData &subset, const std::vector<int> subset_idx);

		std::vector<int> index() { return index_; }
		const std::vector<int>& index() const { return index_; }

		int index(int i) const
		{
			int ni = index_.size();
			if (ni > 0 && i >= ni || i < 0 || i >= number())
				THROW("Aquisition number is out of range");
			if (ni > 0)
				return index_[i];
			else
				return i;
		}

    	/*! 
    		\brief Reader for ISMRMRD::Acquisition from ISMRMRD file. 
      		*	filename_ismrmrd_with_ext:	filename of ISMRMRD rawdata file with .h5 extension.
      		* 
      		* In case the ISMRMRD::Dataset constructor throws an std::runtime_error the reader catches it, 
      		* displays the message and throws it again.
			* To avoid reading noise samples and other calibration data, the TO_BE_IGNORED macro is employed
			* to exclude potentially incompatible input. 
    	*/
		void read( const std::string& filename_ismrmrd_with_ext );

	protected:
		bool sorted_ = false;
		std::vector<int> index_;
        std::vector<KSpaceSorting> sorting_;
		AcquisitionsInfo acqs_info_;

		static std::string _storage_scheme;
		// new MRAcquisitionData objects will be created from this template
		// using same_acquisitions_container()
		static gadgetron::shared_ptr<MRAcquisitionData> acqs_templ_;

		virtual MRAcquisitionData* clone_impl() const = 0;
		MRAcquisitionData* clone_base() const;

	private:
		void binary_op_(int op, 
			const MRAcquisitionData& a_x, const MRAcquisitionData& a_y,
			complex_float_t a = 0, complex_float_t b = 0);

	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief File implementation of Abstract MR acquisition data container class.

	Acquisitions are stored in HDF5 file.
	*/
	class AcquisitionsFile : public MRAcquisitionData {
	public:
		AcquisitionsFile() { own_file_ = false; }
		AcquisitionsFile
			(std::string filename, bool create_file = false,
			AcquisitionsInfo info = AcquisitionsInfo());
		AcquisitionsFile(AcquisitionsInfo info);
		~AcquisitionsFile();

		// initializes MRAcquisitionData template as AcquisitionsFile
		static void init() {
			static bool initialized = false;
			if (!initialized) {
				acqs_templ_.reset(new AcquisitionsFile());
				_storage_scheme = "file";
				MRAcquisitionData::storage_scheme();
				initialized = true;
			}
		}
		// sets MRAcquisitionData template as AcquisitionsFile
		static void set_as_template()
		{
			init();
			acqs_templ_.reset(new AcquisitionsFile);
			_storage_scheme = "file";
		}

		// implements 'overwriting' of an acquisition file data with new values:
		// in reality, creates new file with new data and deletes the old one
		void take_over_impl(AcquisitionsFile& ac);

		void write_acquisitions_info();

		// implementations of abstract methods

		virtual void empty();
		virtual void take_over(MRAcquisitionData& ad)
		{
			AcquisitionsFile& af = dynamic_cast<AcquisitionsFile&>(ad);
			take_over_impl(af);
		}
		virtual void set_data(const complex_float_t* z, int all = 1);
		virtual unsigned int items() const;
		virtual unsigned int number() const { return items(); }
		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const;
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
		{
			//std::cerr << 
			THROW("AcquisitionsFile::set_acquisition not implemented yet, sorry\n");
		}
		virtual void append_acquisition(ISMRMRD::Acquisition& acq);
		virtual void copy_acquisitions_info(const MRAcquisitionData& ac);
		virtual void copy_acquisitions_data(const MRAcquisitionData& ac);

		virtual AcquisitionsFile*
			same_acquisitions_container(const AcquisitionsInfo& info) const
		{
			return (AcquisitionsFile*) new AcquisitionsFile(info);
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			init();
			DataContainer* ptr = acqs_templ_->same_acquisitions_container(acqs_info_);
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(ptr));
		}
		virtual gadgetron::unique_ptr<MRAcquisitionData> new_acquisitions_container()
		{
			init();
			return gadgetron::unique_ptr<MRAcquisitionData>
				(acqs_templ_->same_acquisitions_container(acqs_info_));
		}

	private:
		bool own_file_;
		std::string filename_;
		gadgetron::shared_ptr<ISMRMRD::Dataset> dataset_;
		virtual AcquisitionsFile* clone_impl() const
		{
			init();
			return (AcquisitionsFile*)clone_base();
		}
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief A vector implementation of the abstract MR acquisition data container
	class.

	Acquisitions are stored in an std::vector<shared_ptr<ISMRMRD::Acquisition> >
	object.
	*/
	class AcquisitionsVector : public MRAcquisitionData {
	public:
		AcquisitionsVector(AcquisitionsInfo info = AcquisitionsInfo())
		{
			acqs_info_ = info;
		}
		static void init() 
		{ 
			AcquisitionsFile::init(); 
		}
		static void set_as_template()
		{
			init();
			acqs_templ_.reset(new AcquisitionsVector);
			_storage_scheme = "memory";
		}
		virtual void empty();
		virtual void take_over(MRAcquisitionData& ad) {}
		virtual unsigned int number() const { return (unsigned int)acqs_.size(); }
		virtual unsigned int items() const { return (unsigned int)acqs_.size(); }
		virtual void append_acquisition(ISMRMRD::Acquisition& acq)
		{
			acqs_.push_back(gadgetron::shared_ptr<ISMRMRD::Acquisition>
				(new ISMRMRD::Acquisition(acq)));
		}
		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const
		{
			int ind = index(num);
			acq = *acqs_[ind];
		}
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
		{
			int ind = index(num);
			*acqs_[ind] = acq;
		}
		virtual void copy_acquisitions_info(const MRAcquisitionData& ac)
		{
			acqs_info_ = ac.acquisitions_info();
		}
		virtual void copy_acquisitions_data(const MRAcquisitionData& ac);
		virtual void set_data(const complex_float_t* z, int all = 1);

		virtual AcquisitionsVector* same_acquisitions_container
			(const AcquisitionsInfo& info) const
		{
			return new AcquisitionsVector(info);
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			init();
			DataContainer* ptr = acqs_templ_->same_acquisitions_container(acqs_info_);
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(ptr));
		}
		virtual gadgetron::unique_ptr<MRAcquisitionData>
			new_acquisitions_container()
		{
			init();
			return gadgetron::unique_ptr<MRAcquisitionData>
				(acqs_templ_->same_acquisitions_container(acqs_info_));
		}

	private:
		std::vector<gadgetron::shared_ptr<ISMRMRD::Acquisition> > acqs_;
		virtual AcquisitionsVector* clone_impl() const
		{
			init();
			return (AcquisitionsVector*)clone_base();
		}
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract Gadgetron image data container class.

	*/

	class ISMRMRDImageData : public MRImageData {
	public:
		//ISMRMRDImageData(ISMRMRDImageData& id, const char* attr, 
		//const char* target); //does not build, have to be in the derived class
		
		virtual void empty() = 0;
		virtual unsigned int number() const = 0;
		virtual gadgetron::shared_ptr<ImageWrap> sptr_image_wrap
			(unsigned int im_num) = 0;
		virtual gadgetron::shared_ptr<const ImageWrap> sptr_image_wrap
			(unsigned int im_num) const = 0;
//		virtual ImageWrap& image_wrap(unsigned int im_num) = 0;
//		virtual const ImageWrap& image_wrap(unsigned int im_num) const = 0;
		virtual void append(int image_data_type, void* ptr_image) = 0;
		virtual void append(const ImageWrap& iw) = 0;
		virtual void get_data(complex_float_t* data) const;
		virtual void set_data(const complex_float_t* data);
		virtual void get_real_data(float* data) const;
		virtual void set_real_data(const float* data);
		virtual int read(std::string filename, std::string variable = "", int iv = -1);
		virtual void write(const std::string &filename, const std::string &groupname, const bool dicom) const;
		virtual void write(const std::string &filename) const { this->write(filename, "", false); }
		virtual Dimensions dimensions() const
		{
			Dimensions dim;
			const ImageWrap& iw = image_wrap(0);
			int d[4];
			iw.get_dim(d);
			dim["x"] = d[0];
			dim["y"] = d[1];
			dim["z"] = d[2];
			dim["c"] = d[3];
			dim["n"] = number();
			return dim;
		}
        virtual void get_image_dimensions(unsigned int im_num, int* dim) const
		{
			if (im_num >= number())
				dim[0] = dim[1] = dim[2] = dim[3] = 0;
            const ImageWrap&  iw = image_wrap(im_num);
			iw.get_dim(dim);
		}
		virtual gadgetron::shared_ptr<ISMRMRDImageData> 
			new_images_container() const = 0;
		virtual gadgetron::shared_ptr<ISMRMRDImageData>
			clone(const char* attr, const char* target) = 0;
		virtual int image_data_type(unsigned int im_num) const
		{
			return image_wrap(im_num).type();
		}

		virtual float norm() const;
		virtual void dot(const DataContainer& dc, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		virtual void multiply(
			const DataContainer& a_x,
			const DataContainer& a_y);
		virtual void divide(
			const DataContainer& a_x,
			const DataContainer& a_y);

		virtual void sort() = 0;
		bool sorted() const { return sorted_; }
		void set_sorted(bool sorted) { sorted_ = sorted; }
		std::vector<int> index() { return index_; }
		const std::vector<int>& index() const { return index_; }
		int index(int i) const
		{
			int ni = index_.size();
			if (ni > 0 && i >= ni || i < 0 || i >= number())
				THROW("Image number is out of range");
			if (ni > 0)
				return index_[i];
			else
				return i;
		}
		ImageWrap& image_wrap(unsigned int im_num)
		{
			gadgetron::shared_ptr<ImageWrap> sptr_iw = sptr_image_wrap(im_num);
			return *sptr_iw;
		}
		const ImageWrap& image_wrap(unsigned int im_num) const
		{
			const gadgetron::shared_ptr<const ImageWrap>& sptr_iw = 
				sptr_image_wrap(im_num);
			return *sptr_iw;
		}
        /// Set the meta data
        void set_meta_data(const AcquisitionsInfo &acqs_info);
        /// Get the meta data
        const AcquisitionsInfo &get_meta_data() const { return acqs_info_; }


	protected:
		bool sorted_=false;
		std::vector<int> index_;
        AcquisitionsInfo acqs_info_;
	};

	typedef ISMRMRDImageData GadgetronImageData;

	/*!
	\ingroup Gadgetron Data Containers
	\brief A vector implementation of the abstract Gadgetron image data 
	container class.

	Images are stored in an std::vector<shared_ptr<ImageWrap> > object.
	*/

	class GadgetronImagesVector : public GadgetronImageData {
	public:
		typedef ImageData::Iterator BaseIter;
		typedef ImageData::Iterator_const BaseIter_const;
		typedef std::vector<gadgetron::shared_ptr<ImageWrap> >::iterator
			ImageWrapIter;
		typedef std::vector<gadgetron::shared_ptr<ImageWrap> >::const_iterator 
			ImageWrapIter_const;
		class Iterator : public MRImageData::Iterator {
		public:
			Iterator(ImageWrapIter iw, int n, int i, const ImageWrap::Iterator& it) :
				iw_(iw), n_(n), i_(i), iter_(it), end_((**iw).end())
			{}
			Iterator(const Iterator& iter) : iw_(iter.iw_), n_(iter.n_), i_(iter.i_),
				iter_(iter.iter_), end_(iter.end_), sptr_iter_(iter.sptr_iter_)
			{}
			Iterator& operator=(const Iterator& iter)
			{
				iw_ = iter.iw_;
				n_ = iter.n_;
				i_ = iter.i_;
				iter_ = iter.iter_;
				end_ = iter.end_;
				sptr_iter_ = iter.sptr_iter_;
				return *this;
			}
			virtual bool operator==(const BaseIter& ai) const
			{
				DYNAMIC_CAST(const Iterator, i, ai);
				return iter_ == i.iter_;
			}
			virtual bool operator!=(const BaseIter& ai) const
			{
				DYNAMIC_CAST(const Iterator, i, ai);
				return iter_ != i.iter_;
			}
			Iterator& operator++()
			{
				if (i_ >= n_ || (i_ == n_ - 1 && iter_ == end_))
					throw std::out_of_range("cannot advance out-of-range iterator");
				++iter_;
				if (iter_ == end_ && i_ < n_ - 1) {
					++i_;
					++iw_;
					iter_ = (**iw_).begin();
					end_ = (**iw_).end();
				}
				return *this;
			}
			Iterator& operator++(int)
			{
				sptr_iter_.reset(new Iterator(*this));
				if (i_ >= n_ || (i_ == n_ - 1 && iter_ == end_))
					throw std::out_of_range("cannot advance out-of-range iterator");
				++iter_;
				if (iter_ == end_ && i_ < n_ - 1) {
					++i_;
					++iw_;
					iter_ = (**iw_).begin();
					end_ = (**iw_).end();
				}
				return *sptr_iter_;
			}
			NumRef& operator*()
			{
				if (i_ >= n_ || (i_ == n_ - 1 && iter_ == end_))
					throw std::out_of_range
					("cannot dereference out-of-range iterator");
				return *iter_;
			}
		private:
			//std::vector<gadgetron::shared_ptr<ImageWrap> >::iterator iw_;
			ImageWrapIter iw_;
			int n_;
			int i_;
			ImageWrap::Iterator iter_;
			ImageWrap::Iterator end_;
			gadgetron::shared_ptr<Iterator> sptr_iter_;
		};

		class Iterator_const : public MRImageData::Iterator_const {
		public:
			Iterator_const(ImageWrapIter_const iw, int n, int i, 
				const ImageWrap::Iterator_const& it) :
				iw_(iw), n_(n), i_(i), iter_(it), end_((**iw).end_const())
			{}
			Iterator_const(const Iterator_const& iter) : iw_(iter.iw_),
				n_(iter.n_), i_(iter.i_),
				iter_(iter.iter_), end_(iter.end_), sptr_iter_(iter.sptr_iter_)
			{}
			Iterator_const& operator=(const Iterator_const& iter)
			{
				iw_ = iter.iw_;
				n_ = iter.n_;
				i_ = iter.i_;
				iter_ = iter.iter_;
				end_ = iter.end_;
				sptr_iter_ = iter.sptr_iter_;
                return *this;
			}
			bool operator==(const BaseIter_const& ai) const
			{
				DYNAMIC_CAST(const Iterator_const, i, ai);
				return iter_ == i.iter_;
			}
			bool operator!=(const BaseIter_const& ai) const
			{
				DYNAMIC_CAST(const Iterator_const, i, ai);
				return iter_ != i.iter_;
			}
			Iterator_const& operator++()
			{
				if (i_ >= n_ || (i_ == n_ - 1 && iter_ == end_))
					throw std::out_of_range("cannot advance out-of-range iterator");
				++iter_;
				if (iter_ == end_ && i_ < n_ - 1) {
					++i_;
					++iw_;
					iter_ = (**iw_).begin_const();
					end_ = (**iw_).end_const();
				}
				return *this;
			}
			// causes crashes and very inefficient anyway
			//Iterator_const& operator++(int)
			//{
			//	sptr_iter_.reset(new Iterator_const(*this));
			//	if (i_ >= n_ || (i_ == n_ - 1 && iter_ == end_))
			//		throw std::out_of_range("cannot advance out-of-range iterator");
			//	++iter_;
			//	if (iter_ == end_ && i_ < n_ - 1) {
			//		++i_;
			//		++iw_;
			//		iter_ = (**iw_).begin_const();
			//		end_ = (**iw_).end_const();
			//	}
			//	return *sptr_iter_;
			//}
			const NumRef& operator*() const
			{
				if (i_ >= n_ || (i_ == n_ - 1 && iter_ == end_))
					throw std::out_of_range
					("cannot dereference out-of-range iterator");
				ref_.copy(*iter_);
				return ref_;
				//return *iter_;
			}
		private:
			//std::vector<gadgetron::shared_ptr<ImageWrap> >::const_iterator iw_;
			ImageWrapIter_const iw_;
			int n_;
			int i_;
			ImageWrap::Iterator_const iter_;
			ImageWrap::Iterator_const end_;
			mutable NumRef ref_;
			gadgetron::shared_ptr<Iterator_const> sptr_iter_;
		};

		GadgetronImagesVector() : images_()
		{}
        GadgetronImagesVector(const GadgetronImagesVector& images);
		GadgetronImagesVector(GadgetronImagesVector& images, const char* attr,
			const char* target);
		virtual void empty()
		{
			images_.clear();
		}
		virtual unsigned int items() const
		{ 
			return (unsigned int)images_.size(); 
		}
		virtual unsigned int number() const 
		{ 
			return (unsigned int)images_.size(); 
		}
		virtual void append(int image_data_type, void* ptr_image)
		{
			images_.push_back(gadgetron::shared_ptr<ImageWrap>
				(new ImageWrap(image_data_type, ptr_image)));
		}
		virtual void append(const ImageWrap& iw)
		{
			images_.push_back(gadgetron::shared_ptr<ImageWrap>(new ImageWrap(iw)));
		}
		virtual void sort();
		virtual gadgetron::shared_ptr<ImageWrap> sptr_image_wrap
			(unsigned int im_num)
		{
			int i = index(im_num);
			return images_.at(i);
		}
		virtual gadgetron::shared_ptr<const ImageWrap> sptr_image_wrap
			(unsigned int im_num) const
		{
			int i = index(im_num);
			return images_.at(i);
		}
/*		virtual ImageWrap& image_wrap(unsigned int im_num)
		{
			gadgetron::shared_ptr<ImageWrap> sptr_iw = sptr_image_wrap(im_num);
			return *sptr_iw;
		}
		virtual const ImageWrap& image_wrap(unsigned int im_num) const
		{
			const gadgetron::shared_ptr<const ImageWrap>& sptr_iw = 
				sptr_image_wrap(im_num);
			return *sptr_iw;
		}
*/
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(new_images_container()));
		}
		virtual gadgetron::shared_ptr<GadgetronImageData> new_images_container() const
		{
			gadgetron::shared_ptr<GadgetronImageData> sptr_img
				((GadgetronImageData*)new GadgetronImagesVector());
			sptr_img->set_meta_data(get_meta_data());
			return sptr_img;
		}
		virtual gadgetron::shared_ptr<GadgetronImageData>
			clone(const char* attr, const char* target)
		{
			return gadgetron::shared_ptr<GadgetronImageData>
				(new GadgetronImagesVector(*this, attr, target));
		}

		virtual Iterator& begin()
		{
			ImageWrapIter iw = images_.begin();
			begin_.reset(new Iterator(iw, images_.size(), 0, (**iw).begin()));
			return *begin_;
		}
		virtual Iterator& end()
		{
			ImageWrapIter iw = images_.begin();
			int n = images_.size();
			for (int i = 0; i < n - 1; i++)
				++iw;
			end_.reset(new Iterator(iw, n, n - 1, (**iw).end()));
			return *end_;
		}
		virtual Iterator_const& begin() const
		{
			ImageWrapIter_const iw = images_.begin();
			begin_const_.reset
				(new Iterator_const(iw, images_.size(), 0, (**iw).begin_const()));
			return *begin_const_;
		}
		virtual Iterator_const& end() const
		{
			ImageWrapIter_const iw = images_.begin();
			int n = images_.size();
			for (int i = 0; i < n - 1; i++)
				++iw;
			end_const_.reset
				(new Iterator_const(iw, n, n - 1, (**iw).end_const()));
			return *end_const_;
		}
		virtual void get_data(complex_float_t* data) const;
		virtual void set_data(const complex_float_t* data);
		virtual void get_real_data(float* data) const;
		virtual void set_real_data(const float* data);

        /// Clone and return as unique pointer.
        std::unique_ptr<GadgetronImagesVector> clone() const
        {
            return std::unique_ptr<GadgetronImagesVector>(this->clone_impl());
        }

        /// Print header info
        void print_header(const unsigned im_num);

        /// Is complex?
        virtual bool is_complex() const;

        /// Reorient image. Requires that dimensions match
        virtual void reorient(const VoxelisedGeometricalInfo3D &geom_info_out);

        /// Populate the geometrical info metadata (from the image's own metadata)
        virtual void set_up_geom_info();

    private:
        /// Clone helper function. Don't use.
        virtual GadgetronImagesVector* clone_impl() const
        {
            return new GadgetronImagesVector(*this);
        }

        std::vector<gadgetron::shared_ptr<ImageWrap> > images_;
        mutable gadgetron::shared_ptr<Iterator> begin_;
        mutable gadgetron::shared_ptr<Iterator> end_;
        mutable gadgetron::shared_ptr<Iterator_const> begin_const_;
        mutable gadgetron::shared_ptr<Iterator_const> end_const_;
    };

    /*!
    \ingroup Gadgetron Data Containers
    \brief A coil sensitivities container based on the GadgetronImagesVector class.

    Coil sensitivities can be computed directly form Acquisition data where in the first
    step the vector of the GadgetronImagesVector is used to store the individual channels.

    Then these images are stored in memory, used to compute the coilmaps themselves, the
    individual channel reconstructions are discareded and the coilmaps are stored in their stead.

    This leaves the possibility to compute coilmaps directly based on readily available indi-
    vidual channel reconstructions.
    */

    class CoilSensitivitiesVector : public GadgetronImagesVector
    {

    public:

        CoilSensitivitiesVector() : GadgetronImagesVector(){}
        CoilSensitivitiesVector(const char * file)
        {
            throw std::runtime_error("This has not been implemented yet.");
        }

        void set_csm_smoothness(int s){csm_smoothness_ = s;}

        void calculate(const MRAcquisitionData& acq)
        {
            this->calculate_images(acq);
            this->calculate_csm();
        }

        void calculate_images(const MRAcquisitionData& acq);
        void calculate_csm(GadgetronImagesVector iv);

        CFImage get_csm_as_cfimage(size_t const i) const;

        void get_dim(size_t const num_csm, int* dim) const
        {
            GadgetronImagesVector::get_image_dimensions(num_csm, dim);

        }

    protected:

        bool flag_imgs_suitable_for_csm_computation_=false;


        void calculate_csm(void);

        void calculate_csm(ISMRMRD::NDArray<complex_float_t>& cm, ISMRMRD::NDArray<float>& img, ISMRMRD::NDArray<complex_float_t>& csm);

        void forward(){
            throw LocalisedException("This has not been implemented yet." , __FILE__, __LINE__);
        }
        void backward(){
            throw LocalisedException("This has not been implemented yet." , __FILE__, __LINE__);
        }


    private:
        int csm_smoothness_=0;
        void smoothen_(int nx, int ny, int nz, int nc, complex_float_t* u, complex_float_t* v, int* obj_mask, int w);
        void mask_noise_(int nx, int ny, int nz, float* u, float noise, int* mask);
        float max_diff_(int nx, int ny, int nz, int nc, float small_grad, complex_float_t* u, complex_float_t* v);
        float max_(int nx, int ny, int nz, float* u);

    };

}

#endif





