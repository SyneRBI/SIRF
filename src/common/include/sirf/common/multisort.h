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

namespace NewMultisort {

	template<typename T>
	bool less(const std::vector<T>& a, const std::vector<T>& b)
	{
		size_t n = std::max(a.size(), b.size());
		for (unsigned int i = 0; i < n; i++) {
			if (a[i] < b[i])
				return true;
			if (a[i] > b[i])
				return false;
		}
		return false; // all equal
	}

	template<typename T>
	void sort(std::vector<std::vector<T> > v, int* index)
	{
		int n = v.size();
		std::iota(index, index + n, 0);
		std::stable_sort
			(index, index + n, [&v](int i, int j){return less(v[i], v[j]); });
	}

} // namespace NewMultisort

#endif