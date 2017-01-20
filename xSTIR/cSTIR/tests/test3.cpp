#include <string>

#include "data_handle.h"
#include "cstir.h"
#include "stir.h"
#include "stir_x.h"
#include "tests.h"

void test3()
{
	TextWriter w;
	openChannel(0, &w);

	std::string path("../../examples/");
	std::string filename;
	int dim[3];
	size_t size, sinos, views, tangs, segments;

	filename = path + "my_forward_projection.hs";
	//filename = "tmp.hs";
	boost::shared_ptr<ProjData> sptr_ad = ProjData::read_from_file(filename);
	boost::shared_ptr<ProjDataInfo> sptr_adi = sptr_ad->get_proj_data_info_sptr();
	size = sptr_ad->size_all();
	segments = sptr_ad->get_num_segments();
	sinos = sptr_ad->get_num_sinograms();
	views = sptr_ad->get_num_views();
	tangs = sptr_ad->get_num_tangential_poss();
	std::cout << "segments: " << segments << '\n';
	std::cout << "sinograms: " << sinos << '\n';
	std::cout << "views: " << views << '\n';
	std::cout << "tangential positions: " << tangs << '\n';
	std::cout << "size: " << size << ' ' << sinos*views*tangs << '\n';

	//const ProjDataInfo& adi = *sptr_adi;
	//Voxels3DF voxels(adi);
	//Voxels3DF voxels(*sptr_adi);
	Voxels3DF* ptr_voxels = new Voxels3DF(*sptr_adi);
	sptrImage3DF sptr_image(ptr_voxels);
	Voxels3DF& voxels = *ptr_voxels;
	Coord3DF vsize = voxels.get_voxel_size();
	std::cout << vsize[1] << '\n';
	std::cout << vsize[2] << '\n';
	std::cout << vsize[3] << '\n';
	//sptrImage3DF sptr_image(voxels.clone());
	//filename = path + "my_image.hv";
	//sptrImage3DF sptr_image(read_from_file<Image3DF>(filename));
	Image3DF& image = *sptr_image;
	image.fill(1.0f);
	xSTIR_getImageDimensions(image, dim);
	std::cout << dim[0] << '\n';
	std::cout << dim[1] << '\n';
	std::cout << dim[2] << '\n';

	OBJECT(ProjMatrixByBin, RayTracingMatrix, matrix, sptr_matrix);
	PETAcquisitionModelUsingMatrix<Image3DF> acq_mod;
	acq_mod.set_matrix(sptr_matrix);
	acq_mod.set_up(sptr_ad, sptr_image);

	//boost::shared_ptr<ProjData> sptr_fd = acq_mod.forward(image);

	//double* acq_data = new double[size];
	//sptr_ad->copy_to(acq_data);
	//double* fwd_data = new double[size];
	//sptr_fd->copy_to(fwd_data);
	//std::cout << "acq diff: " << diff(size, acq_data, fwd_data) << '\n';

	//delete[] acq_data;
	//delete[] fwd_data;
}