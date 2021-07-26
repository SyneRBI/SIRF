#include <cctype>
#include <string>

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
}

