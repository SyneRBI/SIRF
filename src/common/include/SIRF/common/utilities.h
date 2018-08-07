#ifndef SIRF_UTILITIES
#define SIRF_UTILITIES

#include <chrono>
#include <string>

namespace sirf {

	class SIRFUtilities {
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
			long long int ms = milliseconds();
			calls++;
			sprintf(buff, "tmp_%d_%lld.h5", calls, ms);
			return std::string(buff);
		}

	};

}

#endif