#pragma once

#ifndef CASE_INSENSITIVE_STRING_COMPARISON
#define CASE_INSENSITIVE_STRING_COMPARISON

#include <string>

/*!
\file
\ingroup Common
\brief Case insensitive string comparison, replaces boost::iequals.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/

namespace sirf {
	bool iequals(const std::string& a, const std::string& b);
}

#endif