#pragma once

#include <memory>

namespace sirf {
	template<class Vector>
	class Operator {
	public:
		virtual std::shared_ptr<Vector> apply(const Vector& v) = 0;
		std::shared_ptr<Vector> operator()(const Vector& v)
		{
			return apply(v);
		}
	};
}
