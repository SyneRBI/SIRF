#include <cctype>
#include <string>

#include "sirf/common/iequals.h"

namespace sirf {
	bool iequals(const std::string& a, const std::string& b)
	{
		unsigned int n = a.size();
		if (b.size() != n)
			return false;
		for (unsigned int i = 0; i < n; i++)
			if (tolower(a[i]) != tolower(b[i]))
				return false;
		return true;
	}
	std::string path(const std::string& path)
	{
		std::string output = path;
		for (int i = 0; i < path.size(); i++)
			if (path[i] == '\\')
				output[i] = '/';
		return output;
	}
}

