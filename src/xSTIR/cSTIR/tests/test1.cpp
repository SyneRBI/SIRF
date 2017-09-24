#include <string>

#include "data_handle.h"
#include "cstir.h"
#include "stir.h"
#include "stir_x.h"
#include "tests.h"
//#include "xstir.h"

double diff(size_t n, double* u, double*v)
{
	double umax = 0;
	double vmax = 0;
	double d = 0;
	for (size_t i = 0; i < n; i++) {
		double ui = abs(u[i]);
		double vi = abs(v[i]);
		double di = abs(u[i] - v[i]);
		if (ui > umax)
			umax = ui;
		if (vi > vmax)
			vmax = vi;
		if (di > d)
			d = di;
	}
	double uvmax = (umax > vmax ? umax : vmax);
	return d / uvmax;
}

double dot(size_t n, double* u, double*v)
{
	double d = 0;
	for (size_t i = 0; i < n; i++)
		d += u[i] * v[i];
	return d;
}

//int xSTIR_getImageDimensions(const Image3DF& image, int* dim)
//{
//	dim[0] = 0;
//	dim[1] = 0;
//	dim[2] = 0;
//	Coordinate3D<int> min_indices;
//	Coordinate3D<int> max_indices;
//	if (!image.get_regular_range(min_indices, max_indices))
//		return -1;
//	for (int i = 0; i < 3; i++)
//		dim[i] = max_indices[i + 1] - min_indices[i + 1] + 1;
//	return 0;
//}

int xSTIR_getImageDataAsDoubleArray(const Image3DF& image, double* data)
{
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		return -1;
	for (int z = min_indices[1], i = 0; z <= max_indices[1]; z++) {
		for (int y = min_indices[2]; y <= max_indices[2]; y++) {
			for (int x = min_indices[3]; x <= max_indices[3]; x++, i++) {
				data[i] = image[z][y][x];
			}
		}
	}
	return 0;
}

void test1()
{
	try {
		TextWriter w;
		openChannel(0, &w);

		std::string path("../../examples/");
		std::string filename;
		int dim[3];
		size_t size, sinos, views, tangs, segments;

		//filename = path + "expected_image.hv";
		filename = path + "my_image.hv";
		sptrImage3DF sptr_image(read_from_file<Image3DF>(filename));
		Image3DF& image = *sptr_image;
		xSTIR_getImageDimensions(image, dim);
		size_t image_size = dim[0] * dim[1] * dim[2];
		std::cout << dim[0] << '\n';
		std::cout << dim[1] << '\n';
		std::cout << dim[2] << '\n';
		double* image_data = new double[image_size];
		xSTIR_getImageDataAsDoubleArray(image, image_data);

		OBJECT(ProjMatrixByBin, RayTracingMatrix, matrix, sptr_matrix);
		//matrix.set_num_tangential_LORs(2);

		//OBJECT(ProjectorByBinPair, ProjectorPairUsingMatrix, ppm, sptr_ppm);
		//ppm.set_proj_matrix_sptr(sptr_matrix);

		filename = path + "my_forward_projection.hs";
		//filename = "tmp.hs";
		std::shared_ptr<ProjData> sptr_ad = ProjData::read_from_file(filename);
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

		std::shared_ptr<ProjDataInfo> sptr_pdi = sptr_ad->get_proj_data_info_sptr();
		ProjDataInfoCylindrical* ptr_pdic = (ProjDataInfoCylindrical*)sptr_pdi.get();
		double rs = ptr_pdic->get_ring_spacing();
		std::cout << "ring spacing: " << rs << '\n';

		double* acq_data = new double[size];
		sptr_ad->copy_to(acq_data);

		std::shared_ptr<ProjData> sptr_a(
			new ProjDataInMemory(sptr_ad->get_exam_info_sptr(),
			sptr_ad->get_proj_data_info_sptr()));
		sptr_a->fill(0.05f);

		std::shared_ptr<ProjData> sptr_b(
			new ProjDataInMemory(sptr_ad->get_exam_info_sptr(),
			sptr_ad->get_proj_data_info_sptr()));
		sptr_b->fill(0.05f);

		std::shared_ptr<ProjData> sptr_nd(
			new ProjDataInMemory(sptr_ad->get_exam_info_sptr(),
			sptr_ad->get_proj_data_info_sptr()));
		sptr_nd->fill(2.0f);
		std::shared_ptr<BinNormalisation> sptr_n(
			new BinNormalisationFromProjData(sptr_nd));
		sptr_n->set_up(sptr_ad->get_proj_data_info_sptr());

		//PETAcquisitionModel<Image3DF> acq_mod;
		//acq_mod.set_up(sptr_ppm, sptr_ad, sptr_image);
		PETAcquisitionModelUsingMatrix<Image3DF> acq_mod;
		acq_mod.set_matrix(sptr_matrix);
		acq_mod.set_up(sptr_ad, sptr_image);
		std::shared_ptr<ProjectorByBinPair> sptr_ppm = acq_mod.projectors_sptr();

		std::shared_ptr<ProjData> sptr_fd = acq_mod.forward(image, "tmp1.hs");
		size = sptr_fd->size_all();
		segments = sptr_fd->get_num_segments();
		sinos = sptr_fd->get_num_sinograms();
		views = sptr_fd->get_num_views();
		tangs = sptr_fd->get_num_tangential_poss();
		std::cout << "segments: " << segments << '\n';
		std::cout << "sinograms: " << sinos << '\n';
		std::cout << "views: " << views << '\n';
		std::cout << "tangential positions: " << tangs << '\n';
		std::cout << "size: " << size << ' ' << sinos*views*tangs << '\n';

		double* fwd_data = new double[size];
		sptr_fd->copy_to(fwd_data);
		std::cout << "acq diff: " << diff(size, acq_data, fwd_data) << '\n';

		sptrImage3DF sptr_im = acq_mod.backward(*sptr_fd);
		Image3DF& im = *sptr_im;
		xSTIR_getImageDimensions(im, dim);
		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		size_t im_size = dim[0] * dim[1] * dim[2];
		double* im_data = new double[im_size];
		xSTIR_getImageDataAsDoubleArray(im, im_data);

		std::cout << dot(size, acq_data, fwd_data) << '\n';
		std::cout << dot(im_size, image_data, im_data) << '\n';

		delete[] acq_data;
		delete[] fwd_data;
		delete[] im_data;

		acq_mod.set_additive_term(sptr_a);
		//acq_mod.set_background_term(sptr_b);
		//acq_mod.set_normalisation(sptr_n);
		sptr_nd->fill(0.5f);
		acq_mod.set_normalisation(sptr_nd);
		sptr_nd->fill(2.0f);
		sptr_fd = acq_mod.forward(image);

		OBJECT(Prior3DF, QuadPrior3DF, prior, sptr_prior);
		prior.set_penalisation_factor(0.00f);

		OBJECT(DataProcessor3DF, CylindricFilter3DF, filter, sptr_filter);

		OBJECT(ObjectiveFunction3DF, PoissonLogLhLinModMeanProjData3DF, obj_fun,
			sptr_fun);
		obj_fun.set_zero_seg0_end_planes(true);
		//obj_fun.set_max_segment_num_to_process(4);
		//obj_fun.set_max_segment_num_to_process(3);
		obj_fun.set_projector_pair_sptr(sptr_ppm);
		obj_fun.set_proj_data_sptr(sptr_fd);
		obj_fun.set_prior_sptr(sptr_prior);
		obj_fun.set_additive_proj_data_sptr(sptr_a);
		obj_fun.set_normalisation_sptr(sptr_n);

		int num_subiterations = 2;
		OBJECT(Reconstruction<Image3DF>, OSMAPOSLReconstruction<Image3DF>, recon,
		//OBJECT(xSTIR_IterativeReconstruction3DF, 
		//	OSMAPOSLReconstruction<Image3DF>, recon,
			sptr_recon);
		recon.set_MAP_model("multiplicative");
		recon.set_num_subsets(12);
		recon.set_num_subiterations(num_subiterations);
		recon.set_save_interval(num_subiterations);
		recon.set_inter_iteration_filter_interval(1);
		recon.set_output_filename_prefix("reconstructedImage");
		recon.set_objective_function_sptr(sptr_fun);
		recon.set_inter_iteration_filter_ptr(sptr_filter);

		//Succeeded s = 
		//	//recon.set_up(sptr_image);
		//	xSTIR_setupReconstruction((void*)&sptr_recon, sptr_image);
		//if (s != Succeeded::yes) {
		//	std::cout << "xSTIR_setupReconstruction failed\n";
		//}
		xSTIR_IterativeReconstruction3DF& xrecon =
			(xSTIR_IterativeReconstruction3DF&)(recon);
		//(xSTIR_IterativeReconstruction3DF&)(*sptr_recon);
		Succeeded s = Succeeded::no;
		if (!xrecon.post_process()) {
			s = xrecon.setup(sptr_image);
			xrecon.subiteration() = xrecon.get_start_subiteration_num();
		}

		for (int iter = 0; iter < num_subiterations; iter++) {
			std::cout << "iteration " << iter << '\n';
			xrecon.update(image);
			//xSTIR_updateReconstruction((void*)&sptr_recon, image);
		}

		std::cout << dim[0] << ' ' << dim[1] << ' ' << dim[2] << '\n';
		double* rec_image_data = new double[image_size];
		xSTIR_getImageDataAsDoubleArray(image, rec_image_data);
		std::cout << "images diff: " << diff(image_size, image_data, rec_image_data) << '\n';
		delete[] rec_image_data;
		delete[] image_data;
	}
	catch (...) {
		std::cout << "exception thrown\n";
	}
}