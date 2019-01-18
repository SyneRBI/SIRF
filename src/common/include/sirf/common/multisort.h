#ifndef MULTISORT
#define MULTISORT

#include <algorithm> // stable_sort
#include <array> // array
#include <numeric> // iota
#include <vector> // vector

namespace Multisort {

	template<typename T, size_t N>
	bool less(const std::array<T, N>& a, const std::array<T, N>& b)
	{
		for (unsigned int i = 0; i < N; i++) {
			if (a[i] < b[i])
				return true;
			if (a[i] > b[i])
				return false;
		}
		return false; // all equal
	}

	template<typename T, size_t N>
	void sort(std::vector<std::array<T, N> > v, int* index)
	{
		int n = v.size();
		std::iota(index, index + n, 0);
		std::stable_sort
			(index, index + n, [&v](int i, int j){return less(v[i], v[j]); });
	}

} // namespace Multisort

#endif