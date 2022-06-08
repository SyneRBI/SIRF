#pragma once

#include <memory>

namespace sirf {
	template<class Vector>
	class Operator {
	public:
		virtual std::shared_ptr<Vector> apply(Vector& v) = 0;
		std::shared_ptr<Vector> operator()(Vector& v)
		{
			return apply(v);
		}
	};
}