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
		float rightmost(Operator<vector_type>& A, vector_type& x)
		{
			value_type lmd;
			value_type a[] = { 0, 0, 0, 0 };
			value_type mu[2];
			value_type u[2];
			value_type v[2];
			value_type zero = 0;
			value_type t;

			std::unique_ptr<vector_type> uptr_y = x.clone();
			vector_type& y = *uptr_y;
			std::unique_ptr<vector_type> uptr_z = x.clone();
			vector_type& z = *uptr_z;
			std::unique_ptr<vector_type> uptr_w = x.clone();
			vector_type& w = *uptr_w;
			std::unique_ptr<vector_type> uptr_Az = x.clone();
			vector_type& Az = *uptr_Az;

			float s = sqrt(abs(x.dot(x)));
			x.scale(s);
			std::shared_ptr<vector_type> sptr_Ax = A(x);
			vector_type& Ax = *sptr_Ax;

			for (int it = 0; it < nit_; it++) {
				lmd = Ax.dot(x);
				y.axpby(1.0, Ax, -lmd, x);
				//std::cout << it << ": " << real_(lmd) << '\n';
				if (it) {
					t = Az.dot(z);
					w.axpby(1.0, Az, -t, z);
					t = w.dot(y) / (lmd - t);
					y.axpby(1.0, y, t, z);
				}
				s = sqrt(abs(y.dot(y)));
				//std::cout << s << '\n';
				if (s == 0.0)
					break;
				y.scale(s);
				t = y.dot(x);
				y.axpby(1.0, y, -t, x);
				s = sqrt(abs(y.dot(y)));
				if (s == 0.0)
					break;
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
				s = sqrt(abs(x.dot(x)));
				//std::cout << s << '\n';
				x.scale(s);
			}
			return real_(lmd);
		}
	private:
		int nit_;
		float real_(value_type z) const
		{
			value_type zero = 0;
			value_type two = 2;
			value_type zc = (z == zero ? 0 : abs(z*z) / z);
			value_type zr = (z + zc)/two;
			float r = abs(zr);
			if (abs(z - r) < abs(z + r))
				return r;
			else
				return -r;
		}
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
			value_type pc = (p == zero ? 0 : abs(p*p) / p);
			value_type qc = (q == zero ? 0 : abs(q*q) / q);
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