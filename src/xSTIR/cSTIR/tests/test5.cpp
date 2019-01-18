#include <iostream>

#include "stir/common.h"
#include "stir/IO/stir_ecat_common.h"

#include "object.h"
#include "sirf/cSTIR/stir_x.h"
#include "sirf/cGadgetron/gadgetron_data_containers.h"

using namespace stir;
using namespace ecat;
using namespace sirf;

int test_a(shared_ptr<ProjData> sptr_data, shared_ptr<Image3DF>& sptr_image);
int test_b(const PETAcquisitionData& acq_data, STIRImageData& image);
int test_c(const ImageData& image);

int test5()
{
	std::string filename;
	int status;
	int dim[4];
	size_t sinos, views, tangs;
	float im_norm;

	std::string SIRF_path = std::getenv("SIRF_PATH");
	if (SIRF_path.length() < 1) {
		std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
		return 1;
	}

	filename = SIRF_path + "/data/examples/PET/my_forward_projection.hs";

	CREATE_OBJECT(PETAcquisitionData, PETAcquisitionDataInFile,
		acq_data, sptr_ad, filename.c_str());
	sinos = acq_data.get_num_sinograms();
	views = acq_data.get_num_views();
	tangs = acq_data.get_num_tangential_poss();
	float acq_norm = acq_data.norm();
	std::cout << "sinograms: " << sinos << '\n';
	std::cout << "views: " << views << '\n';
	std::cout << "tangential positions: " << tangs << '\n';
	std::cout << "acquisition data norm: " << acq_norm << '\n';
	shared_ptr<ProjData> sptr_data = acq_data.data();
	shared_ptr<Image3DF> sptr_image;

	status = test_a(sptr_data, sptr_image);
	if (status)
		return status;
	STIRImageData image(sptr_image);
	image.get_dimensions(dim);
	std::cout << "image dimensions: "
		<< dim[0] << 'x' << dim[1] << 'x' << dim[2] << '\n';
	im_norm = image.norm();
	std::cout << "image norm: " << im_norm << '\n';

	status = test_b(acq_data, image);
	if (status)
		return status;
	image.get_dimensions(dim);
	std::cout << "image dimensions: "
		<< dim[0] << 'x' << dim[1] << 'x' << dim[2] << '\n';
	im_norm = image.norm();
	std::cout << "image norm: " << im_norm << '\n';

	std::cout << "\ntesting conversion from PET data...\n";
	status = test_c(image);
	if (status)
		return status;

	std::cout << "\ntesting conversion from MR data...\n";
	GadgetronImagesVector mr_image;
	mr_image.read(SIRF_path + "/examples/Python/MR/Gadgetron/output.h5");
	mr_image.get_image_dimensions(0, dim);
	int ni = mr_image.number();
	std::cout << ni << " MR images of dimensions nx = "
		<< dim[0] << ", ny = " << dim[1] << ", nz = " << dim[2] << ", nc = " << dim[3] << '\n';
	status = test_c(mr_image);
	if (status)
		return status;

	std::cout << "Press any key to continue";
	getc(stdin);
	return status;
}

// STIR test
int test_a(shared_ptr<ProjData> sptr_data, shared_ptr<Image3DF>& sptr_image)
{
	try {
		FBP2DReconstruction fbp2d;
		fbp2d.set_input_data(sptr_data);
		sptr_image.reset(fbp2d.construct_target_image_ptr());
		fbp2d.reconstruct(sptr_image);
	}
	catch (...) {
		std::cout << "exception thrown\n";
		return 1;
	}
	return 0;
}

// SIRF test
int test_b(const PETAcquisitionData& acq_data, STIRImageData& image)
{
	try {
		xSTIR_FBP2DReconstruction fbp2d;

		fbp2d.set_input(acq_data);
		fbp2d.process();
		shared_ptr<STIRImageData> sptr_image = fbp2d.get_output();

		image.set_data_sptr(sptr_image->data_sptr());
	}
	catch (...) {
		std::cout << "exception thrown\n";
		return 1;
	}
	return 0;
}

int test_c(const ImageData& image)
{
	try {
		STIRImageData img(image);
		int dim[3];
		float im_norm;
		img.get_dimensions(dim);
		std::cout << "image dimensions: "
			<< dim[0] << 'x' << dim[1] << 'x' << dim[2] << '\n';
		im_norm = img.norm();
		std::cout << "image norm: " << im_norm << '\n';
	}
	catch (...) {
		std::cout << "exception thrown\n";
		return 1;
	}
	return 0;
}

int test_d()
{
	try {
	}
	catch (...) {
		std::cout << "exception thrown\n";
		return 1;
	}
	return 0;
}