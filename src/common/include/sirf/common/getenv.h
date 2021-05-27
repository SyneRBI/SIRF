#pragma once

#ifndef SIRF_GETENV
#define SIRF_GETENV

#include "sirf/iUtilities/LocalisedException.h"

namespace sirf {
	std::string getenv(const char* name)
	{
		const char* value = std::getenv(name);
		std::string s;
		if (value)
			s = value;
		else
			THROW(s + "??? Environmental variable " + name + " not defined\n");
		return s;
	}
}

#endif