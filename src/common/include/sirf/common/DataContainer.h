/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC

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

#ifndef SIRF_ABSTRACT_DATA_CONTAINER_TYPE
#define SIRF_ABSTRACT_DATA_CONTAINER_TYPE

#include<complex>
#include <map>
#include "sirf/iUtilities/DataHandle.h"

/*!
\ingroup Common
\brief Abstract data container.

A class for a set of generally heterogeneous items of data.

Has vector features: norm, dot product, linear combination,
which rely on the same features of the items.
*/

namespace sirf {

	typedef std::map<std::string, int> Dimensions;

	class DataContainer {
	public:
		virtual ~DataContainer() {}
		//virtual DataContainer* new_data_container() const = 0;
		virtual ObjectHandle<DataContainer>* new_data_container_handle() const = 0;
		virtual unsigned int items() const = 0;
		virtual bool is_complex() const = 0;
		/// returns the size of data elements
		virtual int bits() const
		{
			// default value
			return is_complex() ? 16 * sizeof(float) : 8 * sizeof(float);
		}

		/// returns the norm of this container viewed as a vector
		virtual float norm() const = 0;

		/// below all void* are actually either float* (STIR containers and NiftiImageData)
		/// or complex_float_t* (Gadgetron containers)

		/// calculates the dot product of this container with another one
		virtual void dot(const DataContainer& dc, void* ptr) const = 0;

		/// calculates the sum of this container elements
		virtual void sum(void* ptr) const = 0;

		/// calculates the value of this container's element with the largest real part
		virtual void max(void* ptr) const = 0;

		/// calculates the value of this container's element with the smallest real part
		virtual void min(void* ptr) const = 0;

		/// \c *this = the elementwise product \c x*y
		virtual void multiply
			(const DataContainer& x, const DataContainer& y) = 0;
		/// \c *this = the product \c x * y with scalar y
		virtual void multiply
			(const DataContainer& x, const void* ptr_y) = 0;

		/// \c *this = the sum \c x + y with scalar y
		virtual void add
			(const DataContainer& x, const void* ptr_y) = 0;

		/// \c *this = the elementwise ratio \c x / y
		virtual void divide
			(const DataContainer& x, const DataContainer& y) = 0;

		/// \c *this = the elementwise \c max(x, y)
		virtual void maximum
			(const DataContainer& x, const DataContainer& y) = 0;
		virtual void maximum
			(const DataContainer& x, const void* ptr_y) = 0;

		/// \c *this = the elementwise \c min(x, y)
		virtual void minimum
			(const DataContainer& x, const DataContainer& y) = 0;
		virtual void minimum
			(const DataContainer& x, const void* ptr_y) = 0;

		/// \c *this = the elementwise \c pow(x, y)
		virtual void power
			(const DataContainer& x, const DataContainer& y) = 0;
		virtual void power
			(const DataContainer& x, const void* ptr_y) = 0;

		/// \c *this = the elementwise \c exp(x)
		virtual void exp(const DataContainer& x) = 0;
		/// \c *this = the elementwise \c log(x)
		virtual void log(const DataContainer& x) = 0;
		/// \c *this = the elementwise \c sqrt(x)
		virtual void sqrt(const DataContainer& x) = 0;
		/// \c *this = the elementwise \c sign(x)
		virtual void sign(const DataContainer& x) = 0;
		/// \c *this = the elementwise \c abs(x)
		virtual void abs(const DataContainer& x) = 0;

		/// \c *this = the linear combination of \c x and \c y
		virtual void axpby(
			const void* ptr_a, const DataContainer& x,
			const void* ptr_b, const DataContainer& y) = 0;
		/// alternative interface to the above
		virtual void xapyb(
			const DataContainer& x, const void* ptr_a,
			const DataContainer& y, const void* ptr_b) = 0;

		/// \c *this = elementwise sum of two elementwise products \c x*a and \c y*b
		virtual void xapyb(
			const DataContainer& x, const DataContainer& a,
			const DataContainer& y, const DataContainer& b) = 0;

		/// \c *this = elementwise sum of \c x*a and elementwise \c y*b
		virtual void xapyb(
			const DataContainer& a_x, const void* ptr_a,
			const DataContainer& a_y, const DataContainer& a_b) = 0;

		/// \c *this = elementwise sum of elementwise \c x*a and \c y*b
		void xapyb(
			const DataContainer& a_x, const DataContainer& a_a,
			const DataContainer& a_y, const void* ptr_b)
		{
			xapyb(a_y, ptr_b, a_x, a_a);
		}

		virtual void write(const std::string &filename) const = 0;

		bool is_empty() const
		{
			return items() < 1;
		}

		std::unique_ptr<DataContainer> clone() const
		{
			return std::unique_ptr<DataContainer>(this->clone_impl());
		}

		/// overwrites this container's complex data with complex conjugate values
		void conjugate()
		{
			this->conjugate_impl();
		}

		///  returns unique pointer to the complex-conjugated copy of this container
		std::unique_ptr<DataContainer> conjugate() const
		{
			DataContainer* ptr = this->clone_impl();
			ptr->conjugate();
			return std::unique_ptr<DataContainer>(ptr);
		}

		template<typename T>
		static T product(T x, T y)
		{
			return x * y;
		}

		template<typename T>
		static T ratio(T x, T y)
		{
			return x / y;
		}

		template<typename T>
		static T inverse_ratio(T x, T y)
		{
			return y / x;
		}

		template<typename T>
		static T sum(T x, T y)
		{
			return x + y;
		}

		template<typename T>
		static T maximum(T x, T y)
		{
			return std::max(x, y);
		}
		template<typename T>
		static T maxabs(T x, T y)
		{
			return std::max(std::abs(x), std::abs(y));
		}
		template<typename T>
		static T maxreal(T x, T y)
		{
			return std::real(x) > std::real(y) ? x : y;
		}

		template<typename T>
		static T minimum(T x, T y)
		{
			return std::min(x, y);
		}
		template<typename T>
		static T minabs(T x, T y)
		{
			return std::min(std::abs(x), std::abs(y));
		}
		template<typename T>
		static T minreal(T x, T y)
		{
			return std::real(x) < std::real(y) ? x : y;
		}
		static std::complex<float> power(std::complex<float> x, std::complex<float> y)
		{
			return std::pow(x, y);
		}
		static std::complex<float> exp(std::complex<float> x)
		{
			return std::exp(x);
		}
		static std::complex<float> log(std::complex<float> x)
		{
			return std::log(x);
		}
		static std::complex<float> sqrt(std::complex<float> x)
		{
			return std::sqrt(x);
		}
		template<typename T>
		static T sign(T x)
		{
			return (std::real(x) > 0) - (std::real(x) < 0);
		}
		template<typename T>
		static T abs(T x)
		{
			return T(std::abs(x));
		}

	protected:
		virtual DataContainer* clone_impl() const = 0;
		/// we assume data to be real, complex data containers must override this
		virtual void conjugate_impl()
		{
			if (is_complex())
				THROW("complex data containes must override conjugate_impl()");
		}
	};
}

#endif
