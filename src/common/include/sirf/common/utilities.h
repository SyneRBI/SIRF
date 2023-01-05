#ifndef SIRF_UTILITIES
#define SIRF_UTILITIES

namespace sirf {
	char path_separator();
    std::string append_path(const std::string& path, int na, const char* a, ...);
    std::string append_path(std::string path, const char* a, ...);
}

#endif