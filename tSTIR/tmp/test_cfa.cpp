#include <math.h>
#include <stdio.h>

void compute_flat_areas
(int nx, int ny, const double* image, double noise, int* flat_area);

double f_image(double sigma, double x, double y)
{
	return exp(-(x*x + y*y) / (sigma*sigma));
}

void print_array(int nx, int ny, int array[])
{
	for (int iy = 0; iy < ny; iy++) {
		for (int ix = 0; ix < nx; ix++)
			printf("%3d", array[ix + nx*iy]);
		printf("\n");
	}
}

void test_cfa()
{
	const int nx = 20;
	const int ny = 20;
	int ix, iy;
	double hx, hy;
	double sigma;
	double noise;
	double image[nx*ny];
	int fa[nx*ny];

	hx = 1.0 / nx;
	hy = 1.0 / ny;
	sigma = 0.4;
	noise = 0.05;

	for (iy = 0; iy < ny; iy++)
		for (ix = 0; ix < nx; ix++)
			image[ix + nx*iy] = f_image(sigma, ix*hx, iy*hy);

	//for (;;) {
	//	printf("x, y: ");
	//	scanf_s("%d %d", &ix, &iy);
	//	if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
	//		break;
	//	printf("%f\n", image[ix + nx*iy]);
	//}

	compute_flat_areas(nx, ny, image, noise, fa);

	//for (;;) {
	//	printf("x, y: ");
	//	scanf_s("%d %d", &ix, &iy);
	//	if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
	//		break;
	//	printf("%d\n", fa[ix + nx*iy]);
	//}
	print_array(nx, ny, fa);
}
