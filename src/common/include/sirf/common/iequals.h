#pragma once

#ifndef CASE_INSENSITIVE_STRING_COMPARISON
#define CASE_INSENSITIVE_STRING_COMPARISON

#include <string>

namespace sirf {
	bool iequals(const std::string& a, const std::string& b);
}

#endif