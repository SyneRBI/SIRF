#ifndef EXTRA_STIR_TYPES
#define EXTRA_STIR_TYPES

#include "cstir.h"
#include "stir.h"
#include "stir_p.h"

#define CAST_PTR(T, X, Y) T* X = (T*)Y

class xSTIR_IterativeReconstruction3DF :
	public IterativeReconstruction<Image3DF> {
public:
	bool post_process() {
		return post_processing();
	}
	Succeeded setup(sptrImage3DF const& image) {
		return set_up(image);
	}
	void update(Image3DF &image) {
		update_estimate(image);
		end_of_iteration_processing(image);
		subiteration_num++;
	}
	int& subiteration() {
		return subiteration_num;
	}
	int subiteration() const {
		return subiteration_num;
	}
	void set_initial_estimate_file(const char* filename) {
		initial_data_filename = filename;
	}
};

class xSTIR_OSSPSReconstruction3DF : public OSSPSReconstruction < Image3DF > {
public:
	float& relaxation_parameter_value() {
		return relaxation_parameter;
	}
};

class xSTIR_GeneralisedObjectiveFunction3DF :
	public GeneralisedObjectiveFunction<Image3DF> {
public:
	bool post_process() {
		return post_processing();
	}
};

class xSTIR_PoissonLogLikelihoodWithLinearModelForMeanAndProjData3DF :
	public PoissonLogLikelihoodWithLinearModelForMeanAndProjData<Image3DF> {
public:
	void set_input_file(const char* filename) {
		input_filename = filename;
	}
};

class xSTIR_GeneralisedPrior3DF : public GeneralisedPrior < Image3DF > {
public:
	bool post_process() {
		return post_processing();
	}
};

class xSTIR_QuadraticPrior3DF : public QuadraticPrior < float > {
public:
	void only2D(int only) {
		only_2D = only != 0;
	}
};

//class xSTIR_ProjectorByBinPairUsingProjMatrixByBin : 
//	public ProjectorByBinPairUsingProjMatrixByBin {
//public:
//	bool post_process() {
//// inaccessible
//		return post_processing();
//	}
//};

inline bool xSTIR_setupPrior(void* ptr)
{
	CAST_PTR(boost::shared_ptr< GeneralisedPrior<Image3DF> >, sptr_obj, ptr);
	CAST_PTR(xSTIR_GeneralisedPrior3DF, obj, sptr_obj->get());
	bool status = obj->post_process();
	return status;
}

inline bool xSTIR_setupObjectiveFunction(void* ptr)
{
	CAST_PTR(boost::shared_ptr< GeneralisedObjectiveFunction<Image3DF> >,
		sptr_obj, ptr);
	CAST_PTR(xSTIR_GeneralisedObjectiveFunction3DF, obj, sptr_obj->get());
	bool status = obj->post_process();
	return status;
}

inline Succeeded xSTIR_setupReconstruction(void* ptr, sptrImage3DF const& image)
{
	CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);
	// not needed - default is non-zero string ("1") anyway
	//recon->set_initial_estimate_file("dummy.hv");
	Succeeded s = Succeeded::no;
	if (recon->post_process())
		return s;
	s = recon->setup(image);
	recon->subiteration() = recon->get_start_subiteration_num();
	return s;
}

inline void xSTIR_updateReconstruction(void* ptr, Image3DF& image) 
{
	CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);
	recon->update(image);
}

inline int& xSTIR_subiteration(void* ptr) 
{
	CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);
	return recon->subiteration();
}

inline void xSTIR_set_initial_estimate_file(void* ptr, const char* filename) 
{
	CAST_PTR(xSTIR_IterativeReconstruction3DF, recon, ptr);
	recon->set_initial_estimate_file(filename);
}

#endif