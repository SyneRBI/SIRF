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
	}
	catch (...) {
		std::cout << "exception thrown\n";
	}
}