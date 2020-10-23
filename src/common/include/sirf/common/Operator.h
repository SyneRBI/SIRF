#pragma once

namespace sirf {
	template<class Vector>
	class Operator {
	public:
		virtual std::shared_ptr<Vector> apply(const Vector& v) const = 0;
		std::shared_ptr<Vector> operator()(const Vector& v) const
		{
			return apply(v);
		}
	};
}