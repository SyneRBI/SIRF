#pragma once

#ifndef DEPRECATION_MACROS
#define DEPRECATION_MACROS

// Deprecation function. With C++14, could use [[deprecated("some message")]]
#if defined(__GNUC__) || defined(__clang__)
#define SIRF_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define SIRF_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define SIRF_DEPRECATED
#endif
#if defined(_MSC_VER)
#define SIRF_DEPRECATED_USING
#else
#define SIRF_DEPRECATED_USING SIRF_DEPRECATED
#endif

#endif