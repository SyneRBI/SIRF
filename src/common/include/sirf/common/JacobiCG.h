#pragma once

#include <cmath>
#include <iostream>

#include <boost/shared_ptr.hpp>

#include "sirf/common/Operator.h"

namespace sirf {
	template<class value_type>
	class JacobiCG {
	public:
		JacobiCG() : nit_(10) {}
		void set_num_iterations(int nit)
		{
			nit_ = nit;
		}
		template<class vector_type>
		value_type rightmost(const Operator<vector_type>& A, vector_type& x) const
		{
			value_type lmd;
			value_type a[] = { 0, 0, 0, 0 };
			value_type mu[2];
			value_type u[2];
			value_type v[2];

			std::unique_ptr<vector_type> uptr_y = x.clone();
			vector_type& y = *uptr_y;
			std::unique_ptr<vector_type> uptr_z = x.clone();
			vector_type& z = *uptr_z;
			std::unique_ptr<vector_type> uptr_w = x.clone();
			vector_type& w = *uptr_w;
			std::unique_ptr<vector_type> uptr_Az = x.clone();
			vector_type& Az = *uptr_Az;

			value_type s = sqrt(x.dot(x));
			x.scale(s);
			std::shared_ptr<vector_type> sptr_Ax = A(x);
			vector_type& Ax = *sptr_Ax;

			for (int it = 0; it < nit_; it++) {
				lmd = Ax.dot(x);
				y.axpby(1.0, Ax, -lmd, x);
				//std::cout << it << ": " << lmd << '\n';
				if (it) {
					s = Az.dot(z);
					w.axpby(1.0, Az, -s, z);
					s = w.dot(y) / (lmd - s);
					y.axpby(1.0, y, s, z);
				}
				s = sqrt(y.dot(y));
				if (s == 0.0)
					return lmd;
				y.scale(s);
				s = y.dot(x);
				y.axpby(1.0, y, -s, x);
				s = sqrt(y.dot(y));
				if (s == 0.0)
					return lmd;
				y.scale(s);
				std::shared_ptr<vector_type> sptr_Ay = A(y);
				vector_type& Ay = *sptr_Ay;
				a[0] = lmd;
				a[1] = Ay.dot(x);
				a[2] = x.dot(Ay);
				a[3] = Ay.dot(y);
				eigh2_(a, mu, u, v);
				z.axpby(u[0], x, u[1], y);
				x.axpby(v[0], x, v[1], y);
				Az.axpby(u[0], Ax, u[1], Ay);
				Ax.axpby(v[0], Ax, v[1], Ay);
				lmd = mu[1];
			}
			return lmd;
		}
	private:
		int nit_;
		void eigh2_(const value_type a[4], value_type lmd[2], value_type x[2], value_type y[2]) const
		{
			value_type zero = 0;
			value_type one = 1;
			value_type c = abs(a[1]);
			value_type d = abs(a[0] - a[3]);
			d = d * d + (c + c) * (c + c);
			value_type s = sqrt(d);
			value_type t = a[0] + a[3];
			if (s == zero && t == zero) {
				lmd[0] = 0.0;
				lmd[1] = 0.0;
				x[0] = 1.0;
				x[1] = 0.0;
				x[2] = 0.0;
				x[3] = 1.0;
				return;
			}
			lmd[0] = (t - s) / (one + one);
			lmd[1] = (t + s) / (one + one);
			value_type p = a[1];
			value_type q = lmd[0] - a[0];
			value_type pc = p == zero ? 0 : abs(p*p) / p;
			value_type qc = q == zero ? 0 : abs(q*q) / q;
			x[0] = p;
			x[1] = q;
			s = sqrt(abs(pc*x[0] + qc*x[1]));
			x[0] /= s;
			x[1] /= s;
			y[0] = -qc / s;
			y[1] = pc / s;
			return;
		}
	};
}