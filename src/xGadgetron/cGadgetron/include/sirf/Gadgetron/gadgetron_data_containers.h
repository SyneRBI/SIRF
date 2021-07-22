/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2020 University College London
Copyright 2020 - 2021 Physikalisch-Technische Bundesanstalt (PTB)

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
\ingroup MR
\brief Specification file for data container classes for Gadgetron data.

\author Evgueni Ovtchinnikov
\author Johannes Mayer
\author SyneRBI
*/

#ifndef GADGETRON_DATA_CONTAINERS
#define GADGETRON_DATA_CONTAINERS

#include <string>
#include <vector>

//#include <boost/algorithm/string.hpp>

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/dataset.h>

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
\ingroup MR
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
\ingroup MR
\brief Serialized ISMRMRD acquisition header (cf. ismrmrd.h).

*/

namespace sirf {

    class FourierEncoding;
    class CartesianFourierEncoding;
#if GADGETRON_TOOLBOXES_AVAILABLE
    class RPEFourierEncoding;
#endif

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
            have_header_ = true;
		}
		std::string data_;
        mutable ISMRMRD::IsmrmrdHeader header_;
        mutable bool have_header_;
	};

    /*!
    \ingroup Gadgetron Data Containers
    \brief Class to keep track of order in k-space
    *
    * The entirety of data consists of all acquisitions in the container. However,
    * the individual acquisitions belong to different subsets of k-space. These
    * each have a different slice, contrast, repetition etc.
    * This class is used to keep track of what acquisitions belong to which subset.
    *

    */
    class KSpaceSubset
    {
        static int const num_kspace_dims_ = 7 + ISMRMRD::ISMRMRD_Constants::ISMRMRD_USER_INTS;

    public:

        typedef std::array<int, num_kspace_dims_> TagType;
        typedef std::vector<int> SetType;

        KSpaceSubset(){
            for(int i=0; i<num_kspace_dims_; ++i)
                this->tag_[i] = -1;
        }

        KSpaceSubset(TagType tag){
            this->tag_ = tag;
            this->idx_set_ = {};
        }

        KSpaceSubset(TagType tag, SetType idx_set){
            this->tag_ = tag;
            this->idx_set_ = idx_set;
        }

        TagType get_tag(void) const {return tag_;}
        SetType get_idx_set(void) const {return idx_set_;}
        void add_idx_to_set(size_t const idx){this->idx_set_.push_back(idx);}

        bool is_first_set() const {
            bool is_first= (tag_[0] == 0);
            if(is_first)
            {
               for(int dim=2; dim<num_kspace_dims_; ++dim)
                   is_first *= (tag_[dim] == 0);
            }
            return is_first;
        }

        static void print_tag(const TagType& tag);
        static void print_acquisition_tag(ISMRMRD::Acquisition acq);

        //! Function to get k-space dimension tag from an ISMRMRD::Acquisition
        /*!
        * This allows to find out which k-space dimension an Acquisition belongs to.
        */
        static TagType get_tag_from_acquisition(ISMRMRD::Acquisition acq);

        //! Function to get k-space dimension tag from an ISMRMRD::Image
        /*!
        * This allows to find out which k-space dimension the image belongs to.
        */
        static TagType get_tag_from_img(const CFImage& img);

    private:
        //! Tag labelling the dimension of k-space
        /*!
        * An int-array of length 7+ISMRMRD_USER_INTS that labels the dimension of k-space.
        * The order is: [average, slice, contrast, phase, repetition, set, segment, user_ (0,...,ISMRMRD_USER_INTS-1)]
        */
        TagType tag_;

        //! Tag labelling the dimension of k-space
        /*!
        * A vector of ints keeping track which acquisitions belong the the
        * dimension labeled by tag_
        */
        SetType idx_set_;
    };

	/*!
	\ingroup MR
	\brief Abstract MR acquisition data container class.

	*/
	class MRAcquisitionData : public DataContainer {
	public:
		// static methods

		// ISMRMRD acquisitions algebra: acquisitions viewed as vectors of 
		// acquisition data
		// y := a x + b y
		static void axpby
			(complex_float_t a, const ISMRMRD::Acquisition& acq_x,
			complex_float_t b, ISMRMRD::Acquisition& acq_y);
		static void xapyb
			(const ISMRMRD::Acquisition& acq_x, complex_float_t a,
			ISMRMRD::Acquisition& acq_y, complex_float_t b);
		static void xapyb
			(const ISMRMRD::Acquisition& acq_x, const ISMRMRD::Acquisition& acq_a,
				ISMRMRD::Acquisition& acq_y, const ISMRMRD::Acquisition& acq_b);

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

        virtual void set_user_floats(float const * const z, int const idx);

		virtual bool is_complex() const
		{
			return true;
		}

		// acquisition data algebra
		virtual void dot(const DataContainer& dc, void* ptr) const;
		virtual void axpby(
			const void* ptr_a, const DataContainer& a_x,
			const void* ptr_b, const DataContainer& a_y);
		virtual void xapyb(
			const DataContainer& a_x, const DataContainer& a_a,
			const DataContainer& a_y, const DataContainer& a_b);
		virtual void xapyb(
			const DataContainer& a_x, const void* ptr_a,
			const DataContainer& a_y, const void* ptr_b);
		//{
		//	axpby(ptr_a, a_x, ptr_b, a_y);
		//}
		virtual void multiply(const DataContainer& x, const DataContainer& y);
		virtual void divide(const DataContainer& x,	const DataContainer& y);
		virtual void maximum(const DataContainer& x, const DataContainer& y)
		{
			THROW("maximum not defined for MRAcquisitionData");
		}
		virtual void minimum(const DataContainer& x, const DataContainer& y)
		{
			THROW("minimum not defined for MRAcquisitionData");
		}
		virtual float norm() const;

		virtual void write(const std::string &filename) const;

		// regular methods

		AcquisitionsInfo acquisitions_info() const { return acqs_info_; }
		void set_acquisitions_info(std::string info) { acqs_info_ = info; }
        void set_acquisitions_info(const AcquisitionsInfo info) { acqs_info_ = info;}

        ISMRMRD::TrajectoryType get_trajectory_type() const;

		gadgetron::unique_ptr<MRAcquisitionData> clone() const
		{
			return gadgetron::unique_ptr<MRAcquisitionData>(this->clone_impl());
		}

		bool undersampled() const;
        int get_acquisitions_dimensions(size_t ptr_dim) const;
        void get_kspace_dimensions(std::vector<size_t>& dims) const;
		uint16_t get_trajectory_dimensions(void) const;
	
		void sort();
		void sort_by_time();
		bool sorted() const { return sorted_; }
		void set_sorted(bool sorted) { sorted_ = sorted; }

        //! Function to get the indices of the acquisitions belonging to different dimensions of k-space
        /*!
        * All acquisitions belong to only one subset in a multi-dimensional k-space.
        * This function returns a vector of sets of indices belonging to the acquisitions of the individual subsets.
        */
        std::vector<KSpaceSubset::SetType > get_kspace_order() const;

        //! Function to get the all KSpaceSubset's of the MRAcquisitionData
        std::vector<KSpaceSubset> get_kspace_sorting() const { return this->sorting_; }

        //! Function to go through Acquisitions and assign them to their k-space dimension
        /*!
        * All acquisitions belong to only one subset in a multi-dimensional k-space. This function goes through
        * all acquisitions in the container, extracts their subset (i.e. which slice contrast etc.) and stores
        * this information s.t. consisten subsets (i.e. all acquisitions belonging to the same slice) can be
        * extracted.
        */
        void organise_kspace();

		virtual std::vector<int> get_flagged_acquisitions_index(const std::vector<ISMRMRD::ISMRMRD_AcquisitionFlags> flags) const;

        virtual void get_subset(MRAcquisitionData& subset, const std::vector<int> subset_idx) const;
        virtual void set_subset(const MRAcquisitionData &subset, const std::vector<int> subset_idx);

		std::vector<int> index() { return index_; }
		const std::vector<int>& index() const { return index_; }

		int index(int i) const
		{
			const std::size_t ni = index_.size();
			if (i < 0 || (ni > 0 && static_cast<std::size_t>(i) >= ni) || static_cast<unsigned>(i) >= number())
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
        std::vector<KSpaceSubset> sorting_;
		AcquisitionsInfo acqs_info_;

		// new MRAcquisitionData objects will be created from this template
		// using same_acquisitions_container()
		static gadgetron::shared_ptr<MRAcquisitionData> acqs_templ_;

		virtual MRAcquisitionData* clone_impl() const = 0;

	private:
		void binary_op_(int op, 
			const MRAcquisitionData& a_x, const MRAcquisitionData& a_y,
			const void* ptr_a = 0, const void* ptr_b = 0);

	};

	/*!
	\ingroup MR
	\brief A vector implementation of the abstract MR acquisition data container
	class.

	Acquisitions are stored in an std::vector<shared_ptr<ISMRMRD::Acquisition> >
	object.
	*/
	class AcquisitionsVector : public MRAcquisitionData {
	public:
        AcquisitionsVector(const std::string& filename_with_ext)
        {
            this->read(filename_with_ext);
        }

        AcquisitionsVector(const AcquisitionsInfo& info = AcquisitionsInfo())
		{
			acqs_info_ = info;
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
			DataContainer* ptr = new AcquisitionsVector(acqs_info_);
			return new ObjectHandle<DataContainer>
				(gadgetron::shared_ptr<DataContainer>(ptr));
		}
		virtual gadgetron::unique_ptr<MRAcquisitionData>
			new_acquisitions_container()
		{
			return gadgetron::unique_ptr<MRAcquisitionData>
				(new AcquisitionsVector(acqs_info_));
		}

	private:
		std::vector<gadgetron::shared_ptr<ISMRMRD::Acquisition> > acqs_;
		virtual AcquisitionsVector* clone_impl() const;
	};

	/*!
	\ingroup MR
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
		virtual void append(gadgetron::shared_ptr<ImageWrap> sptr_iw) = 0;
		virtual gadgetron::shared_ptr<ISMRMRDImageData> abs() const = 0;
		virtual void clear_data()=0;
		virtual void get_data(complex_float_t* data) const;
		virtual void set_data(const complex_float_t* data);
		virtual void get_real_data(float* data) const;
		virtual void set_real_data(const float* data);
		virtual int read(std::string filename, std::string variable = "", int iv = -1);
		virtual void write(const std::string &filename, const std::string &groupname, const bool dicom) const;
		virtual void write(const std::string &filename) const 
		{
			size_t size = filename.size();
			std::string suff = filename.substr(size - 4, 4);
			if (suff == std::string(".dcm")) {
				std::string prefix = filename.substr(0, size - 4);
				this->write(prefix, "", true);
			}
			else {
				std::string fullname = ensure_ext_(filename, "h5");
				this->write(fullname, "", false);
			}
		}
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
        bool check_dimension_consistency() const
        {
            size_t const num_dims = 4;
            std::vector<int> first_img_dims(num_dims), temp_img_dims(num_dims);

            this->get_image_dimensions(0, &first_img_dims[0]);

            bool dims_match = true;
            for(int i=1; i<number(); ++i)
            {
                this->get_image_dimensions(0, &temp_img_dims[0]);
                dims_match *= (first_img_dims == temp_img_dims);
            }
            return dims_match;
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
		virtual void xapyb(
			const DataContainer& a_x, const void* ptr_a,
			const DataContainer& a_y, const void* ptr_b)
		{
			ComplexFloat_ a(*(complex_float_t*)ptr_a);
			ComplexFloat_ b(*(complex_float_t*)ptr_b);
			xapyb_(a_x, a, a_y, b);
		}
		virtual void xapyb(
			const DataContainer& a_x, const void* ptr_a,
			const DataContainer& a_y, const DataContainer& a_b)
		{
			ComplexFloat_ a(*(complex_float_t*)ptr_a);
			DYNAMIC_CAST(const ISMRMRDImageData, b, a_b);
			xapyb_(a_x, a, a_y, b);
		}
		virtual void xapyb(
			const DataContainer& a_x, const DataContainer& a_a,
			const DataContainer& a_y, const void* ptr_b)
		{
			DYNAMIC_CAST(const ISMRMRDImageData, a, a_a);
			ComplexFloat_ b(*(complex_float_t*)ptr_b);
			xapyb_(a_x, a, a_y, b);
		}
		virtual void xapyb(
			const DataContainer& a_x, const DataContainer& a_a,
			const DataContainer& a_y, const DataContainer& a_b)
		{
			DYNAMIC_CAST(const ISMRMRDImageData, a, a_a);
			DYNAMIC_CAST(const ISMRMRDImageData, b, a_b);
			xapyb_(a_x, a, a_y, b);
		}
		virtual void multiply(const DataContainer& x, const DataContainer& y);
		virtual void divide(const DataContainer& x, const DataContainer& y);
		virtual void maximum(const DataContainer& x, const DataContainer& y)
		{
			THROW("maximum not defined for ISMRMRDImageData");
		}
		virtual void minimum(const DataContainer& x, const DataContainer& y)
		{
			THROW("minimum not defined for ISMRMRDImageData");
		}

		void fill(float s);
		void scale(float s);
		complex_float_t dot(const DataContainer& a_x)
		{
			complex_float_t z;
			dot(a_x, &z);
			return z;
		}
		void axpby(
			complex_float_t a, const DataContainer& a_x,
			complex_float_t b, const DataContainer& a_y)
		{
			axpby(&a, a_x, &b, a_y);
		}
		void xapyb(
			const DataContainer& a_x, complex_float_t a,
			const DataContainer& a_y, complex_float_t b)
		{
			xapyb(a_x, &a, a_y, &b);
		}			
		gadgetron::unique_ptr<ISMRMRDImageData> clone() const
		{
			return gadgetron::unique_ptr<ISMRMRDImageData>(this->clone_impl());
		}

		virtual void sort() = 0;
		bool sorted() const { return sorted_; }
		void set_sorted(bool sorted) { sorted_ = sorted; }
		std::vector<int> index() { return index_; }
		const std::vector<int>& index() const { return index_; }
		int index(int i) const
		{
			const std::size_t ni = index_.size();
			if (i < 0 || (ni > 0 && static_cast<std::size_t>(i) >= ni) || static_cast<unsigned>(i) >= number())
				THROW("Image number is out of range. You tried to look up an image number that is not inside the container.");
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
		/// Clone helper function. Don't use.
		virtual ISMRMRDImageData* clone_impl() const = 0;
		std::string ensure_ext_(std::string name, const char* def_ext) const
		{
			auto found = name.find_last_of("/\\");
			if (found == std::string::npos)
				found = name.find_last_of(".");
			else
				found = name.substr(found + 1).find_last_of(".");
			if (found != std::string::npos)
				return name;
			return name + '.' + def_ext;
		}

	private:
		class ComplexFloat_ {
		public:
			ComplexFloat_(complex_float_t v) : v_(v) {}
			unsigned int number() const
			{
				return 0;
			}
			complex_float_t image_wrap(unsigned int i)
			{
				return v_;
			}
		private:
			complex_float_t v_;
		};

		template<class A, class B>
		void xapyb_(const DataContainer& a_x, A& a, const DataContainer& a_y, B& b)
		{
			DYNAMIC_CAST(const ISMRMRDImageData, x, a_x);
			DYNAMIC_CAST(const ISMRMRDImageData, y, a_y);
			unsigned int nx = x.number();
			unsigned int na = a.number();
			unsigned int ny = y.number();
			unsigned int nb = b.number();
			//std::cout << nx << ' ' << ny << '\n';
			if (nx != ny)
				THROW("ImageData sizes mismatch in axpby");
			if (na > 0 && na != nx)
				THROW("ImageData sizes mismatch in axpby");
			if (nb > 0 && nb != nx)
				THROW("ImageData sizes mismatch in axpby");
			unsigned int n = number();
			if (n > 0) {
				if (n != nx)
					THROW("ImageData sizes mismatch in axpby");
				for (unsigned int i = 0; i < nx; i++)
					image_wrap(i).xapyb(x.image_wrap(i), a.image_wrap(i), 
						y.image_wrap(i), b.image_wrap(i));
			}
			else {
				for (unsigned int i = 0; i < nx; i++) {
					const ImageWrap& u = x.image_wrap(i);
					const ImageWrap& v = y.image_wrap(i);
					ImageWrap w(u);
					w.xapyb(u, a.image_wrap(i), v, b.image_wrap(i));
					append(w);
				}
			}
			this->set_meta_data(x.get_meta_data());
		}

	};

	typedef ISMRMRDImageData GadgetronImageData;

	/*!
	\ingroup MR
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

		/*!
		\ingroup MR
		\brief Constructor for images from MR Acquisition data.

		The images are generated with the dimensions given in the recon-space of the 
		MRAcquisitionData's ISMRMRD header information.
		The geometry information and header of the individual images are populated
		based on all consistent subsets of acquisitions in the MRAcquisition object.
		The images can also be created coil-resolved.
		*/
		GadgetronImagesVector(const MRAcquisitionData& ad, const bool coil_resolved=false);
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
        virtual void append(CFImage& img)
        {
            void* vptr_img = new CFImage(img);
            this->append(7, vptr_img);
        }
		virtual void append(const ImageWrap& iw)
		{
			images_.push_back(gadgetron::shared_ptr<ImageWrap>(new ImageWrap(iw)));
		}
		virtual void append(gadgetron::shared_ptr<ImageWrap> sptr_iw)
		{
			images_.push_back(sptr_iw);
		}
		virtual gadgetron::shared_ptr<GadgetronImageData> abs() const;
		virtual void clear_data()
        {
            std::vector<gadgetron::shared_ptr<ImageWrap> > empty_data;
            images_.swap(empty_data);
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
    \ingroup MR
    \brief A coil images container based on the GadgetronImagesVector class.
    */

    class CoilImagesVector : public GadgetronImagesVector
    {
    public:
        CoilImagesVector() : GadgetronImagesVector(){}
        void calculate(const MRAcquisitionData& ad);
    protected:
		std::unique_ptr<MRAcquisitionData> extract_calibration_data(const MRAcquisitionData& ad) const;
	    gadgetron::shared_ptr<FourierEncoding> sptr_enc_;
    };

    /*!
    \ingroup MR
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
		CoilSensitivitiesVector(MRAcquisitionData& ad) :
			GadgetronImagesVector(ad, true){}
        CoilSensitivitiesVector(const char * file)
        {
            throw std::runtime_error("This has not been implemented yet.");
        }

        void set_csm_smoothness(int s){csm_smoothness_ = s;}

        void calculate(CoilImagesVector& iv);
        void calculate(const MRAcquisitionData& acq)
        {
            CoilImagesVector ci;
			ci.calculate(acq);
            calculate(ci);
        }

        CFImage get_csm_as_cfimage(size_t const i) const;
        CFImage get_csm_as_cfimage(const KSpaceSubset::TagType tag, const int offset) const;


        void get_dim(size_t const num_csm, int* dim) const
        {
            GadgetronImagesVector::get_image_dimensions(num_csm, dim);
        }

        void forward(GadgetronImageData& img, GadgetronImageData& combined_img)const;
        void backward(GadgetronImageData& combined_img, const GadgetronImageData& img)const;

    protected:

        void coilchannels_from_combined_image(GadgetronImageData& img, GadgetronImageData& combined_img) const;
        void combine_images_with_coilmaps(GadgetronImageData& combined_img, const GadgetronImageData& img) const;

        void calculate_csm(ISMRMRD::NDArray<complex_float_t>& cm, ISMRMRD::NDArray<float>& img, ISMRMRD::NDArray<complex_float_t>& csm);

    private:
        int csm_smoothness_ = 0;
        void smoothen_(int nx, int ny, int nz, int nc, complex_float_t* u, complex_float_t* v, int* obj_mask, int w);
        void mask_noise_(int nx, int ny, int nz, float* u, float noise, int* mask);
        float max_diff_(int nx, int ny, int nz, int nc, float small_grad, complex_float_t* u, complex_float_t* v);
        float max_(int nx, int ny, int nz, float* u);

    };

void match_img_header_to_acquisition(CFImage& img, const ISMRMRD::Acquisition& acq);

}

#endif
