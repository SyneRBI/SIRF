#include <iostream>

extern "C"
void
find_edges(int nx, int ny, size_t ptr_u, size_t ptr_w)
{
	int mx = 2 * nx - 1;
	int my = 2 * ny - 1;
	int ng = mx*my;
	double* u = (double*)ptr_u;
	float* w = (float*)ptr_w;
	memset(w, 0, ng*sizeof(float));
	float max_grad = 0.0;
	for (int iy = 0; iy < my; iy++) {
		int jy = iy / 2;
		for (int ix = 0; ix < mx; ix++) {
			int jx = ix / 2;
			int i = ix + mx*iy;
			int j = jx + nx*jy;
			if (ix % 2 == 0 && iy % 2 == 0)
				continue;
			if (ix % 2 == 1 && iy % 2 == 1) {
				float s = std::abs(u[j + nx + 1] - u[j]);
				float t = std::abs(u[j + nx] - u[j + 1]);
				if (s > t)
					t = s;
				w[i] = t;
			}
			else if (ix % 2 == 1)
				w[i] = std::abs(u[j + 1] - u[j]);
			else
				w[i] = std::abs(u[j + nx] - u[j]);
			if (w[i] > max_grad)
				max_grad = w[i];
		}
	}
	int nh = ng / 10;
	float* hg = new float[nh];
	memset(hg, 0, nh*sizeof(float));
	for (int iy = 0, i = 0; iy < my; iy++) {
		for (int ix = 0; ix < mx; ix++, i++) {
			int j = (nh - 1)*w[i] / max_grad;
			hg[j]++;
		}
	}
	float cutoff = 0;
	float step = max_grad / nh;
	for (int i = 0; i < nh; i++, cutoff += step)
		if (i > 10 && hg[i] < 1)
			break;

	//std::cout << cutoff << ' ' << max_grad << '\n';

	for (int iy = 0, i = 0; iy < my; iy++) {
		for (int ix = 0; ix < mx; ix++, i++) {
			float t = w[i] / cutoff;
			if (t > 1)
				w[i] = 0.0;
			else
				w[i] = 1.0;
		}
	}
	for (int jy = 0; jy < ny; jy++)
		for (int jx = 0; jx < nx; jx++)
			w[2 * jx + mx * 2 * jy] = 1.0;

	delete[] hg;
}

extern "C"
void
smoothen(int nx, int ny, size_t ptr_u, size_t ptr_w)
{
	int mx = 2 * nx - 1;
	int my = 2 * ny - 1;
	double* u = (double*)ptr_u;
	double* v = new double[nx*ny]; //(double*)ptr_v;
	float* w = (float*)ptr_w;
	for (int jy = 0, j = 0; jy < ny; jy++) {
		int iy = 2 * jy;
		for (int jx = 0; jx < nx; jx++, j++) {
			int ix = 2 * jx;
			int n = 0;
			double r = 0.0;
			double s = 0.0;
			double t = 0.0;
			for (int ky = -1; ky <= 1; ky++) {
				int ly = jy + ky;
				if (ly < 0 || ly >= ny)
					continue;
				for (int kx = -1; kx <= 1; kx++) {
					int lx = jx + kx;
					if (lx < 0 || lx >= nx)
						continue;
					if (kx != 0 || ky != 0) {
						n++;
						r += u[lx + nx*ly];
						s += u[lx + nx*ly] * w[ix + kx + mx*(iy + ky)];
						t += w[ix + kx + mx*(iy + ky)];
					}
				}
			}
			if (t > 0.0)
				v[j] = (u[j] + s / t) / 2;
			else
				v[j] = (u[j] + r / n) / 2;
		}
	}
	memcpy(u, v, nx*ny*sizeof(double));
	delete[] v;
}