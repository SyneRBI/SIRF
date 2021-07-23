/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2020 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

#pragma once

#include <cmath>
#include <iostream>

#include "sirf/common/Operator.h"

/* A simplified implementation of Jacobi-Conjugated Gradients method
   for computing extreme eigenvalues of real symmetric and Hermitian
   operators by E. Ovtchinnikov [SIAM Numer. Anal. 46:2567-2619, 2008].
*/

namespace sirf {
	template<class value_type>
	class JacobiCG {
	public:
		/* This simplified version does fixed number of iterations. */
		JacobiCG() : nit_(10) {}
		void set_num_iterations(int nit)
		{
			nit_ = nit;
		}

		/* Computes the largest eigenvalue of a positive semi-definite operator */
		template<class vector_type>
		float largest(Operator<vector_type>& A,
			vector_type& x /*non-zero initial guess*/, int verb=0)
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
				y.axpby(1.0, Ax, -lmd, x); // residual
				if (verb > 0)
					std::cout << it << ": " << abs(lmd) << '\n';
				if (it) { // conjugate y to the previous search direction
					t = Az.dot(z);
					w.axpby(1.0, Az, -t, z);
					t = w.dot(y) / (lmd - t);
					y.axpby(1.0, y, t, z);
				}
				// normalize y
				s = sqrt(abs(y.dot(y)));
				if (s == 0.0)
					break; // converged
				y.scale(s);
				// orthogonalize y to x
				t = y.dot(x);
				y.axpby(1.0, y, -t, x);
				// normalize y again
				s = sqrt(abs(y.dot(y)));
				if (s == 0.0)
					break; // converged
				y.scale(s);
				// perform Rayleigh-Ritz procedure in span{x,y}
				std::shared_ptr<vector_type> sptr_Ay = A(y);
				vector_type& Ay = *sptr_Ay;
				a[0] = lmd;
				a[1] = Ay.dot(x);
				a[2] = x.dot(Ay);
				a[3] = Ay.dot(y);
				// compute eigenvalues and eigenvectors of 2x2 matrix a
				eigh2_(a, mu, u, v);
				z.axpby(u[0], x, u[1], y);
				x.axpby(v[0], x, v[1], y);
				// span{x,y} = span{x,z} => the new x is a linear combination
				// of the old x and z => we use z as previous search direction
				// on the next iteration
				Az.axpby(u[0], Ax, u[1], Ay);
				Ax.axpby(v[0], Ax, v[1], Ay);
				lmd = mu[1];
//				s = sqrt(abs(x.dot(x)));
				complex_float_t t = x.dot(x);
				s = sqrt(abs(t));
				x.scale(s);
				Ax.scale(s);
			}
			return abs(lmd);
		}
	private:
		int nit_;
		/* computes eigenvalues and eigenvectors of a 2x2 real symmetric
		   or Hermitian matrix [ [a[0], a[1]], [a[2], a[3]] ]
		*/
		void eigh2_(const value_type a[4], value_type lmd[2], value_type x[2],
			value_type y[2]) const
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
