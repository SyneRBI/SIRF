#include <string>

//#include "stir/TextWriter.h"

#include "cstir.h"
#include "data_handle.h"
#include "iutilities.h"
#include "stir.h"
#include "tests.h"

int execution_status(void* handle)
{
	int s = executionStatus(handle);
	if (s)
		std::cout << executionError(handle) << '\n';
	return s;
}

void test2()
{
	std::string path("../../examples/");
	std::string filename;
	int status;
	int dim[3];
	void* handle = 0;
	void* matrix = 0;
	void* image = 0;
	void* img = 0;
	void* am = 0;
	void* ad = 0;
	void* fd = 0;
	void* at = 0;
	void* bt = 0;
	void* nd = 0;
	void* norm = 0;

	TextWriter w;
	openChannel(0, &w);

	for (;;) {
		image = cSTIR_objectFromFile("Image", (path + "my_image.hv").c_str());
		cSTIR_getImageDimensions(image, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

		matrix = cSTIR_newObject("RayTracingMatrix");
		//handle = cSTIR_setParameter
		//	(matrix, "RayTracingMatrix", "num_tangential_LORs", intDataHandle(2));
		//status = execution_status(handle);
		//if (status)
		//	break;
		ad = cSTIR_objectFromFile
			("AcquisitionData", (path + "my_forward_projection.hs").c_str());
		status = execution_status(ad);
		if (status)
			break;
		cSTIR_getAcquisitionsDimensions(ad, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		size_t size = dim[0] * dim[1] * dim[2];

		at = cSTIR_acquisitionsDataFromTemplate(ad);
		cSTIR_getAcquisitionsDimensions(at, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		cSTIR_fillAcquisitionsData(at, 0.05);
		bt = cSTIR_acquisitionsDataFromTemplate(ad);
		cSTIR_getAcquisitionsDimensions(bt, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		cSTIR_fillAcquisitionsData(bt, 0.05);
		nd = cSTIR_acquisitionsDataFromTemplate(ad);
		cSTIR_getAcquisitionsDimensions(nd, (size_t)&dim[0]);
		//std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		cSTIR_fillAcquisitionsData(nd, 2.0);

		am = cSTIR_newObject("PETAcquisitionModelUsingMatrix");
		handle = cSTIR_setParameter(am, "AcqModUsingMatrix", "matrix", matrix);
		status = execution_status(handle);
		if (status)
			break;
		handle = cSTIR_setupAcquisitionModel(am, ad, image);
		status = execution_status(handle);
		if (status)
			break;
		fd = cSTIR_acquisitionModelFwd(am, image, "");
		cSTIR_getAcquisitionsDimensions(fd, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';

		double* adv = new double[size];
		double* fdv = new double[size];

		cSTIR_getAcquisitionsData(ad, (size_t)adv);
		cSTIR_getAcquisitionsData(fd, (size_t)fdv);
		std::cout << "acq diff: " << diff(size, adv, fdv) << '\n';

		img = cSTIR_acquisitionModelBwd(am, fd);
		cSTIR_getImageDimensions(img, (size_t)&dim[0]);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		size_t im_size = dim[0] * dim[1] * dim[2];
		double* image_data = new double[im_size];
		double* img_data = new double[im_size];
		cSTIR_getImageData(image, (size_t)image_data);
		cSTIR_getImageData(img, (size_t)img_data);

		std::cout << dot(size, adv, fdv) << '\n';
		std::cout << dot(im_size, image_data, img_data) << '\n';

		delete[] adv;
		delete[] fdv;

		break;
	}

	deleteDataHandle(handle);
	deleteDataHandle(image);
	deleteDataHandle(matrix);
	deleteDataHandle(am);
	deleteDataHandle(ad);
	deleteDataHandle(fd);
	deleteDataHandle(at);
	deleteDataHandle(bt);
	deleteDataHandle(nd);
}

