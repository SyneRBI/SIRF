#pragma once

#ifndef SIRF_GETENV
#define SIRF_GETENV

namespace sirf {
	std::string getenv(const char* name)
	{
		const char* value = std::getenv(name);
		std::string s;
		if (value)
			s = value;
		return s;
	}
}

#endif