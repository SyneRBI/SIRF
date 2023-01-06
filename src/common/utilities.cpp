#include <cstdarg>
#include <iostream>
#include <sstream>

#include "sirf/common/version.h"

namespace sirf {
	char path_separator()
	{
		std::string filename = __FILE__;
		auto found = filename.find_last_of("/\\");
		if (found == std::string::npos)
			return '/';
		return filename[found];
	}
	std::string append_path(const std::string& path, int na, const char* a, ...)
	{
		va_list args;
		va_start(args, na);
		std::string r = path;
		for (int i = 0; i < na; i++) {
			const char* s = va_arg(args, const char*);
			r += path_separator();
			r += s;
		}
		va_end(args);
		return r;
	}
	std::string append_path(std::string path, const char* a, ...)
	{
		std::string r = path;
		va_list args;
		va_start(args, a);
		r += path_separator();
		r += a;
		for (;;) {
			const char* s = va_arg(args, const char*);
			if (s == NULL)
				break;
			r += path_separator();
			r += s;
		}
		va_end(args);
		return r;
	}
	std::string getenv(const char* name, bool throws = false);
	std::string examples_data_path(const char* data_type)
	{
		std::string SIRF_data_path = sirf::getenv("SIRF_DATA_PATH");
		if (SIRF_data_path.length() > 0)
			return append_path(SIRF_data_path, "examples", data_type, (const char*)NULL);
		std::string SIRF_install_path = sirf::getenv("SIRF_INSTALL_PATH");
		if (SIRF_install_path.length() > 0) {
			std::stringstream sirf_version;
			sirf_version << "SIRF-" << SIRF_VERSION_MAJOR << '.' << SIRF_VERSION_MINOR;
			std::string version = sirf_version.str();
			const char* v = version.c_str();
			return append_path(SIRF_install_path, "share", v, "data", "examples", data_type, (const char*)NULL);
		}
		std::string SIRF_path = sirf::getenv("SIRF_PATH");
		if (SIRF_path.length() > 0)
			return append_path(SIRF_path, "data", "examples", data_type, (const char*)NULL);
		return "";
	}
}