#include <iostream>
#include <cstdarg>

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
}