#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX_AREA 1000

void compute_flat_areas
(int nx, int ny, const double* image, double noise, int* flat_area)
{
	int ix, iy;
	int jx, jy;
	int kx, ky;
	int llist;
	int pass;
	int k, l;
	int* mask;
	int* msk;
	int* x;
	int* y;
	double s, t;

	mask = (int*)malloc(nx*ny*sizeof(int));
	msk = (int*)malloc(nx*ny*sizeof(int));
	x = (int*)malloc(nx*ny*sizeof(int));
	y = (int*)malloc(nx*ny*sizeof(int));

	memset(mask, 0, nx*ny*sizeof(int));
	memset(msk, 0, nx*ny*sizeof(int));

	for (iy = 0; iy < ny; iy++)
		for (ix = 0; ix < nx; ix++) {
			if (mask[ix + nx*iy])
				continue;

			// initialise flat points list
			llist = 1;
			x[0] = ix;
			y[0] = iy;
			msk[ix + nx*iy] = 1;
			mask[ix + nx*iy] = 1;
			s = image[ix + nx*iy];

			// grow flat points list
			for (l = 0; l < llist && llist <= MAX_AREA; l++) {
				kx = x[l];
				ky = y[l];
				for (jy = ky - 1; jy <= ky + 1; jy++) {
					if (jy < 0 || jy >= ny)
						continue;
					for (jx = kx - 1; jx <= kx + 1; jx++) {
						if (jx < 0 || jx >= nx)
							continue;
						//if (mask[jx + nx*jy]) // already listed
						if (msk[jx + nx*jy]) // already listed
							continue;
						t = image[jx + nx*jy];
						if (fabs(s - t) <= noise) {
							msk[jx + nx*jy] = llist;
							mask[jx + nx*jy] = llist;
							x[llist] = jx;
							y[llist] = jy;
							llist++;
						}
					}
				}
			}
			// cleanup
			for (pass = 0; pass < 2; pass++) {
				for (l = 0; l < llist; l++) {
					kx = x[l];
					ky = y[l];
					if (mask[kx + nx*ky] == 0) // removed, skip
						continue;
					// count non-flat neighbours
					k = 0;
					for (jy = ky - 1; jy <= ky + 1; jy++) {
						if (jy < 0 || jy >= ny)
							continue;
						for (jx = kx - 1; jx <= kx + 1; jx++) {
							if (jx < 0 || jx >= nx)
								continue;
							if (msk[jx + nx*jy] == 0) // non-flat neighbour
								k++;
						}
					}
					if (k > 2) // remove from list
						mask[kx + nx*ky] = 0;
				}
			}

			// count remaining points
			k = 0;
			for (l = 0; l < llist; l++) {
				kx = x[l];
				ky = y[l];
				if (mask[kx + nx*ky] == 0) // removed, skip
					continue;
				k++;
			}
			//flat_area[ix + nx*iy] = k;

			for (l = 0; l < llist; l++) {
				kx = x[l];
				ky = y[l];
				flat_area[kx + nx*ky] = k;
				msk[kx + nx*ky] = 0;
			}
		}

	free(mask);
	free(msk);
	free(x);
	free(y);
}

