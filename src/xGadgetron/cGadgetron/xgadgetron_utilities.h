#ifndef XGADGETRON_UTILITIES
#define XGADGETRON_UTILITIES

#include <chrono>
#include <complex>

#include <boost/thread/mutex.hpp>

#include "cgadgetron_shared_ptr.h"

using namespace gadgetron;

class xGadgetronUtilities {
public:
	static long long milliseconds()
	{
		auto now = std::chrono::system_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>
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
	shared_ptr<boost::mutex> sptr()
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
	static shared_ptr<boost::mutex> sptr_mutex_;
	static void init_()
	{
		static bool initialized = false;
		if (!initialized) {
			sptr_mutex_ = shared_ptr<boost::mutex>(new boost::mutex);
			initialized = true;
		}
	}
};

#endif