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
\ingroup xGadgetron Utilities
\brief Various utilities used by SIRF Gadgetron extensions.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef XGADGETRON_UTILITIES
#define XGADGETRON_UTILITIES

#include <chrono>
#include <complex>

#include <boost/thread/mutex.hpp>

#include "sirf/cGadgetron/cgadgetron_shared_ptr.h"

namespace sirf {

	class xGadgetronUtilities {
	public:
		static long long milliseconds()
		{
			auto now = std::chrono::system_clock::now();
			auto ms = std::chrono::duration_cast < std::chrono::milliseconds >
				(now.time_since_epoch());
			return (long long)ms.count();
		}
		static std::string scratch_file_name()
		{
			static int calls = 0;
			char buff[32];
			long long int ms = xGadgetronUtilities::milliseconds();
			calls++;
			sprintf(buff, "tmp_%d_%lld.h5", calls, ms);
			return std::string(buff);
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, unsigned short& t)
		{
			t = (unsigned short)z.real();
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, short& t)
		{
			t = (short)z.real();
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, unsigned int& t)
		{
			t = (unsigned int)z.real();
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, int& t)
		{
			t = (int)z.real();
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, float& t)
		{
			t = (float)z.real();
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, complex_float_t& t)
		{
			t = (complex_float_t)z;
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, double& t)
		{
			t = (double)z.real();
		}
		template<typename T>
		static void convert_complex(std::complex<T> z, complex_double_t& t)
		{
			t = (complex_double_t)z;
		}

	};

	class Mutex {
	public:
		Mutex()
		{
			init_();
		}
		boost::mutex& operator()()
		{
			return *sptr_mutex_.get();
		}
		gadgetron::shared_ptr<boost::mutex> sptr()
		{
			return sptr_mutex_;
		}
		void lock()
		{
			sptr_mutex_->lock();
		}
		void unlock()
		{
			sptr_mutex_->unlock();
		}
	private:
		static gadgetron::shared_ptr<boost::mutex> sptr_mutex_;
		static void init_()
		{
			static bool initialized = false;
			if (!initialized) {
				sptr_mutex_ = gadgetron::shared_ptr<boost::mutex>(new boost::mutex);
				initialized = true;
			}
		}
	};

}

#endif