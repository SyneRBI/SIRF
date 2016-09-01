#include <string>

#include "data_handle.h"
#include "stir.h"
#include "tests.h"

void test1()
{
	try {
		OBJECT(ProjMatrixByBin, RayTracingMatrix, matrix, sptr_matrix);
		matrix.set_num_tangential_LORs(2);

		OBJECT(ProjectorByBinPair, ProjectorPairUsingMatrix, am, sptr_am);
		am.set_proj_matrix_sptr(sptr_matrix);

		std::string path("../../examples/");
		std::string filename(path + "my_forward_projection.hs");
		boost::shared_ptr<ProjData> sptr_ad = ProjData::read_from_file(filename);
		std::cout << "segments: " << sptr_ad->get_num_segments() << '\n';
		size_t size = sptr_ad->size_all();
		size_t sinos = sptr_ad->get_num_sinograms();
		size_t views = sptr_ad->get_num_views();
		size_t tang_pos = sptr_ad->get_num_tangential_poss();
		std::cout << "sinograms: " << sinos << '\n';
		std::cout << "views: " << views << '\n';
		std::cout << "tangential positions: " << tang_pos << '\n';
		std::cout << "size: " << size << ' ' << sinos*views*tang_pos << '\n';
		double* acq_data = new double[size];
		sptr_ad->copy_to(acq_data);
		delete[] acq_data;
	}
	catch (...) {
		std::cout << "exception thrown\n";
	}
}