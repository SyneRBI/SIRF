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

#include <boost/algorithm/string.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>
#include <ismrmrd/meta.h>
#include <ismrmrd/xml.h>

#include "ismrmrd_fftw.h"
#include "cgadgetron_shared_ptr.h"
#include "gadgetron_image_wrap.h"
#include "SIRF/common/data_container.h"
#include "SIRF/common/multisort.h"

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
	class MRAcquisitionData : public aDataContainer < complex_float_t > {
	public:
		MRAcquisitionData() : ordered_(false), index_(0) {}
		virtual ~MRAcquisitionData()
		{
			if (index_)
				delete[] index_;
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
		virtual unsigned int number() = 0;

		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq) = 0;
		virtual void append_acquisition(ISMRMRD::Acquisition& acq) = 0;

		virtual void copy_acquisitions_info(const MRAcquisitionData& ac) = 0;

		// 'export' constructors: workaround for creating 'ABC' objects
		virtual gadgetron::shared_ptr<MRAcquisitionData> new_acquisitions_container() = 0;
		virtual MRAcquisitionData*
			same_acquisitions_container(AcquisitionsInfo info) = 0;

		virtual int set_acquisition_data
			(int na, int nc, int ns, const float* re, const float* im) = 0;

		// acquisition data algebra
		virtual void axpby(
			complex_float_t a, const aDataContainer<complex_float_t>& a_x,
			complex_float_t b, const aDataContainer<complex_float_t>& a_y);
		virtual void multiply(
			const aDataContainer<complex_float_t>& a_x,
			const aDataContainer<complex_float_t>& a_y);
		virtual void divide(
			const aDataContainer<complex_float_t>& a_x,
			const aDataContainer<complex_float_t>& a_y);
		virtual complex_float_t dot(const aDataContainer<complex_float_t>& dc);
		virtual float norm();
		//float diff(MRAcquisitionData& other);

		// regular methods

		std::string acquisitions_info() const { return acqs_info_; }
		void set_acquisitions_info(std::string info) { acqs_info_ = info; }

		bool undersampled() const;
		int get_acquisitions_dimensions(size_t ptr_dim);
		void get_acquisitions_flags(unsigned int n, int* flags);
		unsigned int get_acquisitions_data(unsigned int slice, float* re, float* im);

		void order();
		bool ordered() const { return ordered_; }
		void set_ordered(bool ordered) { ordered_ = ordered; }
		int* index() { return index_; }
		const int* index() const { return index_; }
		int index(int i)
		{
			if (index_ && i >= 0 && i < (int)number())
				return index_[i];
			else
				return i;
		}

		void write(const char* filename);

	protected:
		bool ordered_;
		int* index_;
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

		virtual int set_acquisition_data
			(int na, int nc, int ns, const float* re, const float* im);
		virtual unsigned int items();
		virtual unsigned int number() { return items(); }
		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq);
		virtual void set_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
		{
			std::cerr << "AcquisitionsFile::set_acquisition not implemented yet, sorry\n";
		}
		virtual void append_acquisition(ISMRMRD::Acquisition& acq);
		virtual void copy_acquisitions_info(const MRAcquisitionData& ac);
		virtual MRAcquisitionData*
			same_acquisitions_container(AcquisitionsInfo info)
		{
			return (MRAcquisitionData*) new AcquisitionsFile(info);
		}
		virtual aDataContainer<complex_float_t>*
			new_data_container()
		{
			init();
			return acqs_templ_->same_acquisitions_container(acqs_info_);
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
		virtual unsigned int number() { return (unsigned int)acqs_.size(); }
		virtual unsigned int items() { return (unsigned int)acqs_.size(); }
		virtual void append_acquisition(ISMRMRD::Acquisition& acq)
		{
			acqs_.push_back(gadgetron::shared_ptr<ISMRMRD::Acquisition>
				(new ISMRMRD::Acquisition(acq)));
		}
		virtual void get_acquisition(unsigned int num, ISMRMRD::Acquisition& acq)
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
		virtual int set_acquisition_data
			(int na, int nc, int ns, const float* re, const float* im);
		virtual MRAcquisitionData* same_acquisitions_container(AcquisitionsInfo info)
		{
			return new AcquisitionsVector(info);
		}
		virtual aDataContainer<complex_float_t>* new_data_container()
		{
			AcquisitionsFile::init();
			return acqs_templ_->same_acquisitions_container(acqs_info_);
		}
		virtual gadgetron::shared_ptr<MRAcquisitionData> new_acquisitions_container()
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
	\brief Abstract MR image data container class.

	*/
	class MRImageData : public aDataContainer < complex_float_t > {
	public:
		virtual unsigned int number() = 0;
		virtual int types() = 0;
		virtual void count(int i) = 0;
		virtual gadgetron::shared_ptr<ImageWrap> sptr_image_wrap(unsigned int im_num) = 0;
		virtual gadgetron::shared_ptr<const ImageWrap> sptr_image_wrap
			(unsigned int im_num) const = 0;
		virtual ImageWrap& image_wrap(unsigned int im_num) = 0;
		virtual const ImageWrap& image_wrap(unsigned int im_num) const = 0;
		virtual void append(int image_data_type, void* ptr_image) = 0;
		virtual void append(const ImageWrap& iw) = 0;
		virtual void get_image_dimensions(unsigned int im_num, int* dim) = 0;
		virtual void get_images_data_as_float_array(float* data) const = 0;
		virtual void get_images_data_as_complex_array
			(float* re, float* im) const = 0;
		virtual void set_complex_images_data(const float* re, const float* im) = 0;
		virtual int read(std::string filename) = 0;
		virtual void write(std::string filename, std::string groupname) = 0;
		virtual gadgetron::shared_ptr<MRImageData> new_images_container() = 0;
		virtual gadgetron::shared_ptr<MRImageData>
			clone(unsigned int inc = 1, unsigned int off = 0) = 0;
		virtual gadgetron::shared_ptr<MRImageData>
			clone(const char* attr, const char* target) = 0;
		virtual int image_data_type(unsigned int im_num) const
		{
			return image_wrap(im_num).type();
		}

		virtual void axpby(
			complex_float_t a, const aDataContainer<complex_float_t>& a_x,
			complex_float_t b, const aDataContainer<complex_float_t>& a_y);
		virtual void multiply(
			const aDataContainer<complex_float_t>& a_x,
			const aDataContainer<complex_float_t>& a_y);
		virtual void divide(
			const aDataContainer<complex_float_t>& a_x,
			const aDataContainer<complex_float_t>& a_y);
		virtual complex_float_t dot(const aDataContainer<complex_float_t>& dc);
		virtual float norm();

		void get_image_data_as_cmplx_array
			(unsigned int im_num, float* re, float* im)
		{
			ImageWrap& iw = image_wrap(im_num);
			iw.get_cmplx_data(re, im);
		}

	};

	/*!
	\ingroup Gadgetron Data Containers
	\brief A vector implementation of the abstract MR image data container class.

	Images are stored in an std::vector<shared_ptr<ImageWrap> > object.
	*/
	class ImagesVector : public MRImageData {
	public:
		ImagesVector() : images_(), nimages_(0) {}
		ImagesVector(const ImagesVector& list, const char* attr, const char* target);
		ImagesVector
			(const ImagesVector& list, unsigned int inc = 1, unsigned int off = 0);
		virtual unsigned int items() { return (unsigned int)images_.size(); }
		virtual unsigned int number() { return (unsigned int)images_.size(); }
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
		virtual gadgetron::shared_ptr<ImageWrap> sptr_image_wrap(unsigned int im_num)
		{
			return images_[im_num];
		}
		virtual gadgetron::shared_ptr<const ImageWrap> sptr_image_wrap
			(unsigned int im_num) const
		{
			return images_[im_num];
		}
		virtual ImageWrap& image_wrap(unsigned int im_num)
		{
			gadgetron::shared_ptr<ImageWrap> sptr_iw = sptr_image_wrap(im_num);
			return *sptr_iw;
		}
		virtual const ImageWrap& image_wrap(unsigned int im_num) const
		{
			const gadgetron::shared_ptr<const ImageWrap>& sptr_iw = sptr_image_wrap(im_num);
			return *sptr_iw;
		}
		virtual int read(std::string filename);
		virtual void write(std::string filename, std::string groupname);
		virtual void get_image_dimensions(unsigned int im_num, int* dim)
		{
			if (im_num >= images_.size())
				dim[0] = dim[1] = dim[2] = dim[3] = 0;
			ImageWrap& iw = image_wrap(im_num);
			iw.get_dim(dim);
			//std::string attr = iw.attributes();
			//ISMRMRD::MetaContainer mc;
			//ISMRMRD::deserialize(attr.c_str(), mc);
			//std::cout << mc.as_str("GADGETRON_DataRole") << '\n';
			//std::cout << attr << '\n';
		}
		virtual void get_images_data_as_float_array(float* data) const;
		virtual void get_images_data_as_complex_array(float* re, float* im) const;
		virtual void set_complex_images_data(const float* re, const float* im);
		virtual aDataContainer<complex_float_t>* new_data_container()
		{
			return (aDataContainer<complex_float_t>*)new ImagesVector();
		}
		virtual gadgetron::shared_ptr<MRImageData> new_images_container()
		{
			return gadgetron::shared_ptr<MRImageData>((MRImageData*)new ImagesVector());
		}
		virtual gadgetron::shared_ptr<MRImageData>
			clone(const char* attr, const char* target)
		{
			return gadgetron::shared_ptr<MRImageData>(new ImagesVector(*this, attr, target));
		}
		virtual gadgetron::shared_ptr<MRImageData>
			clone(unsigned int inc = 1, unsigned int off = 0)
		{
			return gadgetron::shared_ptr<MRImageData>(new ImagesVector(*this, inc, off));
		}

	private:
		std::vector<gadgetron::shared_ptr<ImageWrap> > images_;
		int nimages_;

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
	class CoilDataContainer : public aDataContainer < complex_float_t > {
	public:
		virtual float norm()
		{
			return 0.0;
		}
		virtual complex_float_t dot(const aDataContainer<complex_float_t>& dc)
		{
			return complex_float_t(0.0, 0.0);
		}
		virtual void axpby(
			complex_float_t a, const aDataContainer<complex_float_t>& a_x,
			complex_float_t b, const aDataContainer<complex_float_t>& a_y)
		{
			return;
		}
		virtual void multiply(
			const aDataContainer<complex_float_t>& a_x,
			const aDataContainer<complex_float_t>& a_y)
		{
			return;
		}
		virtual void divide(
			const aDataContainer<complex_float_t>& a_x,
			const aDataContainer<complex_float_t>& a_y)
		{
			return;
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
		unsigned int items()
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
		virtual aDataContainer<complex_float_t>* new_data_container()
		{
			return (aDataContainer<complex_float_t>*)new CoilImagesVector();
		}
		virtual unsigned int items()
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

		virtual aDataContainer<complex_float_t>* new_data_container()
		{
			return (aDataContainer<complex_float_t>*)new CoilSensitivitiesAsImages();
		}

		virtual unsigned int items()
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
