#ifndef SIRF_DEPRECATED
#define SIRF_DEPRECATED

namespace sirf {

//! Deprecation macro
#if defined(__GNUC__) || defined(__clang__)
#define SIRF_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define SIRF_DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define SIRF_DEPRECATED
#endif

}

#endif 