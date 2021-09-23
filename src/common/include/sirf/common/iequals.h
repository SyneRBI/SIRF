#pragma once

#ifndef CASE_INSENSITIVE_STRING_COMPARISON
#define CASE_INSENSITIVE_STRING_COMPARISON

#include <string>

/*!
\file
\ingroup Common
\brief Case insensitive string comparison sirf::iequals.

\author Evgueni Ovtchinnikov
\author SyneRBI
*/

namespace sirf {
	/*!
	\ingroup Common
	\brief Case insensitive string comparison, replaces boost::iequals.
	*/
	bool iequals(const std::string& a, const std::string& b);
	void fix_path_separator(std::string& path);
}

#endif