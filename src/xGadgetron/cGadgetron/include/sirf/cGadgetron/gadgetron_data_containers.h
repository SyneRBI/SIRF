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
\ingroup Gadgetron Data Containers
\brief Specification file for data container classes for Gadgetron data.

\author Evgueni Ovtchinnikov
\author CCP PETMR
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
#include "sirf/cGadgetron/ismrmrd_fftw.h"
#include "sirf/cGadgetron/cgadgetron_shared_ptr.h"
#include "sirf/cGadgetron/gadgetron_image_wrap.h"
#include "ismrmrd_fftw.h"

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
	!(acq).isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_SURFACECOILCORRECTIONSCAN_DATA) && \
	(acq).flags() >= (1 << (ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT - 1)))

/*!
\ingroup Gadgetron Data Containers
\brief Serialized ISMRMRD acquisition header (cf. ismrmrd.h).

*/

namespace sirf {

	class AcquisitionsInfo {
	public:
		AcquisitionsInfo(std::string data = "") : data_(data) {}
		AcquisitionsInfo& operator=(std::string data)
		{
			data_ = data;
			return *this;
		}
		const char* c_str() const { return data_.c_str(); }
		operator std::string&() { return data_; }
		operator const std::string&() const { return data_; }

	private:
		std::string data_;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract MR acquisition data container class.

	*/
	class MRAcquisitionData : public DataContainer {
	public:
		MRAcquisitionData() : ordered_(false), index_(0) {}
		virtual ~MRAcquisitionData()
		{
			// if (index_)
				// delete[] index_;
		}

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
		// obsolete
		//static float diff
		//(const ISMRMRD::Acquisition& acq_a, const ISMRMRD::Acquisition& acq_b);

		// abstract methods

		// the number of acquisitions in the container
		virtual unsigned int number() const = 0;

		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const = 0;
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
		virtual void append_acquisition(ISMRMRD::Acquisition& acq) = 0;

		virtual void append_acquisition_sptr( gadgetron::shared_ptr<ISMRMRD::Acquisition> sptr_acqu)=0;
		virtual gadgetron::shared_ptr<ISMRMRD::Acquisition> get_acquisition_sptr(unsigned int num)=0;


		virtual void copy_acquisitions_info(const MRAcquisitionData& ac) = 0;

		// 'export' constructors: workaround for creating 'ABC' objects
		virtual gadgetron::shared_ptr<MRAcquisitionData> new_acquisitions_container() = 0;
		virtual MRAcquisitionData*
			same_acquisitions_container(AcquisitionsInfo info) = 0;

		virtual void set_data(const complex_float_t* z, int all = 1) = 0;

		virtual int set_acquisition_data
			(int na, int nc, int ns, const float* re, const float* im) = 0;

		// acquisition data algebra
		virtual void dot(const DataContainer& dc, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);

		//virtual void axpby(
		//	complex_float_t a, const DataContainer& a_x,
		//	complex_float_t b, const DataContainer& a_y);
		virtual void multiply(
			const DataContainer& a_x,
			const DataContainer& a_y);
		virtual void divide(
			const DataContainer& a_x,
			const DataContainer& a_y);
		//virtual complex_float_t dot(const DataContainer& dc);
		virtual float norm() const;
		//float diff(MRAcquisitionData& other);

		// regular methods

		std::string acquisitions_info() const { return acqs_info_; }
		void set_acquisitions_info(std::string info) { acqs_info_ = info; }

		gadgetron::shared_ptr<MRAcquisitionData> clone();

		bool undersampled() const;
		int get_acquisitions_dimensions(size_t ptr_dim) const;
		//void get_acquisitions_flags(unsigned int n, int* flags);
		unsigned int get_acquisitions_data(unsigned int slice, float* re, float* im) const;
	
		virtual void get_data(complex_float_t* z, int all = 1);

		void order();
		void time_order();
		bool ordered() const { return ordered_; }
		void set_ordered(bool ordered) { ordered_ = ordered; }
		int* index() { return &index_[0]; }
		const int* index() const { return &index_[0]; }
		int index(int i)
		{
			if(i >= 0 && i < (int)number())
			{
				if (index_.size() > 0) 
					return index_[i];
				else
					return i;
			}
			else
			{
				throw LocalisedException("Trying to access out of range acquisition.", __FILE__, __LINE__);
			}
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
		void write(const char* filename);

	protected:
		bool ordered_;
		// int* index_;
		std::vector< int > index_;
		AcquisitionsInfo acqs_info_;

		static std::string _storage_scheme;
		// new MRAcquisitionData objects will be created from this template
		// using same_acquisitions_container()
		static gadgetron::shared_ptr<MRAcquisitionData> acqs_templ_;
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
				initialized = true;
			}
		}
		// sets MRAcquisitionData template as AcquisitionsFile
		static void set_as_template()
		{
			init();
			acqs_templ_.reset(new AcquisitionsFile);
		}

		// implements 'overwriting' of an acquisition file data with new values:
		// in reality, creates new file with new data and deletes the old one
		void take_over(MRAcquisitionData& ac);

		void write_acquisitions_info();

		// implementations of abstract methods

		virtual void set_data(const complex_float_t* z, int all = 1);
		virtual int set_acquisition_data
			(int na, int nc, int ns, const float* re, const float* im);
		virtual unsigned int items() const;
		virtual unsigned int number() const { return items(); }
		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const;
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
		{
			std::cerr << "AcquisitionsFile::set_acquisition not implemented yet, sorry\n";
		}
		virtual void append_acquisition(ISMRMRD::Acquisition& acq);

		virtual void append_acquisition_sptr( gadgetron::shared_ptr<ISMRMRD::Acquisition> sptr_acqu)
		{
			throw LocalisedException("Function undefined for sirf::AcuisitionsFile.", __FILE__, __LINE__);
		};

		virtual gadgetron::shared_ptr<ISMRMRD::Acquisition> get_acquisition_sptr(unsigned int num) 
		{
			throw LocalisedException("Function undefined for sirf::AcuisitionsFile.", __FILE__, __LINE__);
		};

		virtual void copy_acquisitions_info(const MRAcquisitionData& ac);
		virtual MRAcquisitionData*
			same_acquisitions_container(AcquisitionsInfo info)
		{
			return (MRAcquisitionData*) new AcquisitionsFile(info);
		}
		virtual DataContainer*
			new_data_container() const
		{
			init();
			return acqs_templ_->same_acquisitions_container(acqs_info_);
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(new_data_container()));
		}
		virtual gadgetron::shared_ptr<MRAcquisitionData> new_acquisitions_container()
		{
			init();
			return gadgetron::shared_ptr<MRAcquisitionData>
				(acqs_templ_->same_acquisitions_container(acqs_info_));
		}

	private:
		bool own_file_;
		std::string filename_;
		gadgetron::shared_ptr<ISMRMRD::Dataset> dataset_;
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
		static void init() { AcquisitionsFile::init(); }
		static void set_as_template()
		{
			init();
			acqs_templ_.reset(new AcquisitionsVector);
		}
		virtual unsigned int number() const { return (unsigned int)acqs_.size(); }
		virtual unsigned int items() const { return (unsigned int)acqs_.size(); }
		virtual void append_acquisition(ISMRMRD::Acquisition& acq)
		{
			acqs_.push_back(gadgetron::shared_ptr<ISMRMRD::Acquisition>
				(new ISMRMRD::Acquisition(acq)));
		}
		virtual void append_acquisition_sptr( gadgetron::shared_ptr<ISMRMRD::Acquisition> sptr_acqu)
		{
			acqs_.push_back( sptr_acqu );
		}
		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) const
		{
			int ind = index(num);
			acq = *acqs_[ind];
		}
		virtual gadgetron::shared_ptr<ISMRMRD::Acquisition> get_acquisition_sptr(unsigned int num) 
		{
			int ind = index(num);
			return acqs_[ind];
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
		virtual void set_data(const complex_float_t* z, int all = 1);
		virtual int set_acquisition_data
			(int na, int nc, int ns, const float* re, const float* im);
		virtual MRAcquisitionData* same_acquisitions_container
			(AcquisitionsInfo info)
		{
			return new AcquisitionsVector(info);
		}
		virtual DataContainer* new_data_container() const
		{
			AcquisitionsFile::init();
			return acqs_templ_->same_acquisitions_container(acqs_info_);
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(new_data_container()));
		}
		virtual gadgetron::shared_ptr<MRAcquisitionData>
			new_acquisitions_container()
		{
			AcquisitionsFile::init();
			return gadgetron::shared_ptr<MRAcquisitionData>
				(acqs_templ_->same_acquisitions_container(acqs_info_));
		}

	private:
		std::vector<gadgetron::shared_ptr<ISMRMRD::Acquisition> > acqs_;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract Gadgetron image data container class.

	*/

	class ISMRMRDImageData : public MRImageData {
	public:
		ISMRMRDImageData() : ordered_(false), index_(0) {}
		//ISMRMRDImageData(ISMRMRDImageData& id, const char* attr, 
		//const char* target); //does not build, have to be in the derived class
		virtual ~ISMRMRDImageData()
		{
			if (index_)
				delete[] index_;
		}

		virtual unsigned int number() const = 0;
		virtual int types() = 0;
		virtual void count(int i) = 0;
		virtual gadgetron::shared_ptr<ImageWrap> sptr_image_wrap
			(unsigned int im_num) = 0;
		virtual gadgetron::shared_ptr<const ImageWrap> sptr_image_wrap
			(unsigned int im_num) const = 0;
		virtual ImageWrap& image_wrap(unsigned int im_num) = 0;
		virtual const ImageWrap& image_wrap(unsigned int im_num) const = 0;
		virtual void append(int image_data_type, void* ptr_image) = 0;
		virtual void append(const ImageWrap& iw) = 0;
		virtual void get_data(complex_float_t* data) const;
		virtual void set_data(const complex_float_t* data);
		virtual void get_real_data(float* data) const;
		virtual void set_real_data(const float* data);
		virtual int read(std::string filename);
		virtual void write(const std::string &filename, const std::string &groupname) const;
        virtual void write(const std::string &filename) const { this->write(filename,""); }
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
		virtual void get_image_dimensions(unsigned int im_num, int* dim)
		{
			if (im_num >= number())
				dim[0] = dim[1] = dim[2] = dim[3] = 0;
			ImageWrap& iw = image_wrap(im_num);
			iw.get_dim(dim);
		}
		virtual gadgetron::shared_ptr<ISMRMRDImageData> 
			new_images_container() = 0;
		virtual gadgetron::shared_ptr<ISMRMRDImageData>
			clone(const char* attr, const char* target) = 0;
		virtual int image_data_type(unsigned int im_num) const
		{
			return image_wrap(im_num).type();
		}

		virtual void dot(const DataContainer& dc, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		//virtual void axpby(
		//	complex_float_t a, const DataContainer& a_x,
		//	complex_float_t b, const DataContainer& a_y);
		virtual void multiply(
			const DataContainer& a_x,
			const DataContainer& a_y);
		virtual void divide(
			const DataContainer& a_x,
			const DataContainer& a_y);
		//virtual complex_float_t dot(const DataContainer& dc);
		virtual float norm() const;

		void order();
		bool ordered() const { return ordered_; }
		void set_ordered(bool ordered) { ordered_ = ordered; }
		int* index() { return index_; }
		const int* index() const { return index_; }
		int index(int i) const
		{
			if (index_)
				return index_[i];
			else
				return i;
		}

	protected:
		bool ordered_;
		int* index_;
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
				const Iterator& i = (const Iterator&)ai;
				return iter_ == i.iter_;
			}
			virtual bool operator!=(const BaseIter& ai) const
			{
				const Iterator& i = (const Iterator&)ai;
				return iter_ != i.iter_;
			}
			Iterator& operator++()
			{
				if (i_ >= n_ || i_ == n_ - 1 && iter_ == end_)
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
				if (i_ >= n_ || i_ == n_ - 1 && iter_ == end_)
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
				if (i_ >= n_ || i_ == n_ - 1 && iter_ == end_)
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
			}
			bool operator==(const BaseIter_const& ai) const
			{
				const Iterator_const& i = (const Iterator_const&)ai;
				return iter_ == i.iter_;
			}
			bool operator!=(const BaseIter_const& ai) const
			{
				const Iterator_const& i = (const Iterator_const&)ai;
				return iter_ != i.iter_;
			}
			Iterator_const& operator++()
			{
				if (i_ >= n_ || i_ == n_ - 1 && iter_ == end_)
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
			//	if (i_ >= n_ || i_ == n_ - 1 && iter_ == end_)
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
				if (i_ >= n_ || i_ == n_ - 1 && iter_ == end_)
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

		GadgetronImagesVector() : images_(), nimages_(0)
		{}
        GadgetronImagesVector(const GadgetronImagesVector& images);
		GadgetronImagesVector(GadgetronImagesVector& images, const char* attr,
			const char* target);
		virtual unsigned int items() const
		{ 
			return (unsigned int)images_.size(); 
		}
		virtual unsigned int number() const 
		{ 
			return (unsigned int)images_.size(); 
		}
		virtual int types()
		{
			if (nimages_ > 0)
				return (int)(images_.size() / nimages_);
			else
				return 1;
		}
		virtual void count(int i)
		{
			if (i > nimages_)
				nimages_ = i;
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
		virtual gadgetron::shared_ptr<ImageWrap> sptr_image_wrap
			(unsigned int im_num)
		{
			int i = index(im_num);
			return images_[i];
		}
		virtual gadgetron::shared_ptr<const ImageWrap> sptr_image_wrap
			(unsigned int im_num) const
		{
			int i = index(im_num);
			return images_[i];
		}
		virtual ImageWrap& image_wrap(unsigned int im_num)
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
		virtual DataContainer* new_data_container() const
		{
			return (DataContainer*)new GadgetronImagesVector();
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(new_data_container()));
		}
		virtual gadgetron::shared_ptr<GadgetronImageData> new_images_container()
		{
			return gadgetron::shared_ptr<GadgetronImageData>
				((GadgetronImageData*)new GadgetronImagesVector());
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

    protected:
        /// Populate the geometrical info metadata (from the image's own metadata)
        virtual void set_up_geom_info();

	private:
        /// Clone helper function. Don't use.
        virtual GadgetronImagesVector* clone_impl() const
        {
            return new GadgetronImagesVector(*this);
        }

		std::vector<gadgetron::shared_ptr<ImageWrap> > images_;
		int nimages_;
		mutable gadgetron::shared_ptr<Iterator> begin_;
		mutable gadgetron::shared_ptr<Iterator> end_;
		mutable gadgetron::shared_ptr<Iterator_const> begin_const_;
		mutable gadgetron::shared_ptr<Iterator_const> end_const_;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract coil data class.

	Abstract 4-dimensional (x, y, z, coil) single precision complex array.
	*/
	class CoilData {
	public:
		virtual ~CoilData() {}
		virtual void get_dim(int* dim) const = 0;
		virtual void get_data(float* re, float* im) const = 0;
		virtual void set_data(const float* re, const float* im) = 0;
		virtual void get_data(complex_float_t* data) const = 0;
		virtual void set_data(const complex_float_t* data) = 0;
		virtual void get_data_abs(float* v) const = 0;
		virtual complex_float_t& operator()(int x, int y, int z, int c) = 0;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief The ISMRMRD::Image< complex_float_t > implementation
	of the abstract coil data container class.

	*/
	class CoilDataAsCFImage : public CoilData {
	public:
		CoilDataAsCFImage
			(uint16_t nx = 0, uint16_t ny = 1, uint16_t nz = 1, uint16_t nc = 1) :
			img_(nx, ny, nz, nc)
		{
		}
		virtual complex_float_t& operator()(int x, int y, int z, int c)
		{
			return img_(x, y, z, c);
		}
		ISMRMRD::Image < complex_float_t >& image()
		{
			return img_;
		}
		const ISMRMRD::Image < complex_float_t >& image() const
		{
			return img_;
		}
		virtual void get_dim(int* dim) const
		{
			dim[0] = img_.getMatrixSizeX();
			dim[1] = img_.getMatrixSizeY();
			dim[2] = img_.getMatrixSizeZ();
			dim[3] = img_.getNumberOfChannels();
		}
		virtual void get_data(float* re, float* im) const;
		virtual void set_data(const float* re, const float* im);
		virtual void get_data(complex_float_t* data) const
		{
			memcpy(data, img_.getDataPtr(), img_.getDataSize());
		}
		virtual void set_data(const complex_float_t* data)
		{
			memcpy(img_.getDataPtr(), data, img_.getDataSize());
		}
		virtual void get_data_abs(float* v) const;
	private:
		ISMRMRD::Image < complex_float_t > img_;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract coil data container class.

	*/
	class CoilDataContainer : public DataContainer {
	public:
		virtual float norm() const
		{
			return 0.0;
		}
		virtual void dot(const DataContainer& dc, void* ptr) const
		{}
		//virtual complex_float_t dot(const DataContainer& dc)
		//{
		//	return complex_float_t(0.0, 0.0);
		//}
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y)
		{
		}
		//virtual void axpby(
		//	complex_float_t a, const DataContainer& a_x,
		//	complex_float_t b, const DataContainer& a_y)
		//{
		//}
		virtual void multiply(
			const DataContainer& a_x,
			const DataContainer& a_y)
		{
		}
		virtual void divide(
			const DataContainer& a_x,
			const DataContainer& a_y)
		{
		}
		void get_dim(int slice, int* dim) //const
		{
			CoilData& ci = (CoilData&)(*this)(slice);
			ci.get_dim(dim);
		}
		void get_data(int slice, float* re, float* im) //const
		{
			CoilData& ci = (CoilData&)(*this)(slice);
			ci.get_data(re, im);
		}
		void set_data(int slice, float* re, float* im)
		{
			CoilData& ci = (CoilData&)(*this)(slice);
			ci.set_data(re, im);
		}
		void get_data(int slice, complex_float_t* data) //const
		{
			CoilData& ci = (CoilData&)(*this)(slice);
			ci.get_data(data);
		}
		void set_data(int slice, complex_float_t* data)
		{
			CoilData& ci = (CoilData&)(*this)(slice);
			ci.set_data(data);
		}
		void get_data_abs(int slice, float* v) //const
		{
			CoilData& ci = (CoilData&)(*this)(slice);
			ci.get_data_abs(v);
		}
		virtual void append(gadgetron::shared_ptr<CoilData> sptr_csm) = 0;
		virtual CoilData& operator()(int slice) = 0;
		//virtual const CoilData& operator()(int slice) const = 0;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief A vector implementation of the abstract coil data container class.

	CoilData stored in an std::vector<shared_ptr<CoilData> > object.
	*/
	class CoilDataVector {
	public:
		unsigned int items() const
		{
			return (unsigned int)coil_data_.size();
		}
		CoilData& data(int slice) {
			return *coil_data_[slice];
		}
		virtual void append(gadgetron::shared_ptr<CoilData> sptr_cd)
		{
			coil_data_.push_back(sptr_cd);
		}
	private:
		std::vector< gadgetron::shared_ptr<CoilData> > coil_data_;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract coil images container class.

	*/
	class CoilImagesContainer : public CoilDataContainer {
	public:
		virtual CoilData& operator()(int slice) = 0;
		virtual void compute(MRAcquisitionData& ac);
		ISMRMRD::Encoding encoding() const
		{
			return encoding_;
		}
	protected:
		ISMRMRD::Encoding encoding_;
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief A vector implementation of the abstract coil images container class.

	Coil images stored in an std::vector<shared_ptr<CoilData> > object.
	*/
	class CoilImagesVector : public CoilImagesContainer, public CoilDataVector {
	public:
		virtual DataContainer* new_data_container() const
		{
			return (DataContainer*)new CoilImagesVector();
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(new_data_container()));
		}
		virtual unsigned int items() const
		{
			return CoilDataVector::items();
		}
		virtual CoilData& operator()(int slice)
		{
			return data(slice);
		}
		virtual void append(gadgetron::shared_ptr<CoilData> sptr_cd)
		{
			CoilDataVector::append(sptr_cd);
		}
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief Abstract coil sensitivities container class.

	*/
	class CoilSensitivitiesContainer : public CoilDataContainer {
	public:
		void set_csm_smoothness(int s)
		{
			csm_smoothness_ = s;
		}
		virtual CoilData& operator()(int slice) = 0;

		virtual void compute(MRAcquisitionData& ac)
		{
			//if (!ac.ordered())
			//	ac.order();
			CoilImagesVector cis;
			cis.compute(ac);
			compute(cis);
		}

		virtual void compute(CoilImagesContainer& cis);

		void append_csm
			(int nx, int ny, int nz, int nc, const float* re, const float* im)
		{
			//CoilData* ptr_img = new CoilDataType(nx, ny, nz, nc);
			CoilData* ptr_img = new CoilDataAsCFImage(nx, ny, nz, nc);
			gadgetron::shared_ptr<CoilData> sptr_img(ptr_img);
			ptr_img->set_data(re, im);
			append(sptr_img);
		}

	protected:
		int csm_smoothness_;

	private:
		void compute_csm_(
			ISMRMRD::NDArray<complex_float_t>& cm,
			ISMRMRD::NDArray<float>& img,
			ISMRMRD::NDArray<complex_float_t>& csm
			);

		float max_(int nx, int ny, float* u);
		void mask_noise_
			(int nx, int ny, float* u, float noise, int* mask);
		void cleanup_mask_(int nx, int ny, int* mask, int bg, int minsz, int ex);
		void smoothen_
			(int nx, int ny, int nz,
			complex_float_t* u, complex_float_t* v,
			int* obj_mask);
	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief A vector implementation of the abstract coil sensitivities container
	class.

	Coil sensitivities stored in an std::vector<shared_ptr<CoilData> > object.
	*/
	class CoilSensitivitiesAsImages : public CoilSensitivitiesContainer,
		public CoilDataVector {
	public:
		CoilSensitivitiesAsImages()
		{
			csm_smoothness_ = 0;
		}
		CoilSensitivitiesAsImages(const char* file);

		virtual DataContainer* new_data_container() const
		{
			return (DataContainer*)new CoilSensitivitiesAsImages();
		}
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const
		{
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(new_data_container()));
		}

		virtual unsigned int items() const
		{
			return CoilDataVector::items();
		}
		virtual CoilData& operator()(int slice)
		{
			return data(slice);
		}
		virtual void append(gadgetron::shared_ptr<CoilData> sptr_cd)
		{
			CoilDataVector::append(sptr_cd);
		}

	};
}

#endif
