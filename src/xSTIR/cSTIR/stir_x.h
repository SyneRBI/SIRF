/*
CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

#ifndef EXTRA_STIR_TYPES
#define EXTRA_STIR_TYPES

#include <chrono>

#include "data_handle.h"
#include "stir_types.h"

#define MIN_BIN_EFFICIENCY 1.0e-20f

class SIRFUtilities {
public:
	static long long milliseconds()
	{
		auto now = std::chrono::system_clock::now();
		auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
		return (long long)ms.count();
	}
	static std::string scratch_file_name()
	{
		static int calls = 0;
		char buff[32];
		long long int ms = milliseconds();
		calls++;
		sprintf(buff, "tmp_%d_%lld.h5", calls, ms);
		return std::string(buff);
	}
};

template <typename T>
class aDataContainer {
public:
	virtual ~aDataContainer() {}
	virtual boost::shared_ptr<aDataContainer<T> > new_data_container() = 0;
	virtual unsigned int items() = 0;
	virtual float norm() = 0;
	virtual T dot(aDataContainer<T>& dc) = 0;
	virtual void axpby(
		T a, const aDataContainer<T>& x,
		T b, const aDataContainer<T>& y) = 0;
};
#if 0
class PETAcquisitions : public aDataContainer<float> {
public:
	virtual PETAcquisitions* same_acquisitions_container() = 0;
	virtual boost::shared_ptr<PETAcquisitions> new_acquisitions_container() = 0;
protected:
	static boost::shared_ptr<PETAcquisitions> acqs_templ_;
};

class BinAcquisitions : public PETAcquisitions {
public:
	unsigned int items()
	{
		return 1;
	}
	float norm()
	{
		ProjData* self = (ProjData*)this;
		SegmentBySinogram<float> segment = self->get_empty_segment_by_sinogram(0);
		return 0;
	}
	float dot(aDataContainer<float>& other)
	{
		return 0;
	}
	void axpby(float a, const aDataContainer<float>& x, 
		float b, const aDataContainer<float>& y)
	{

	}
};

class BinAcquisitionsInMemory : 
	public BinAcquisitions, public ProjDataInMemory {
public:
	BinAcquisitionsInMemory(shared_ptr<ExamInfo> const& exam_info_sptr,
		shared_ptr<ProjDataInfo> const& proj_data_info_ptr) :
		ProjDataInMemory(exam_info_sptr, proj_data_info_ptr)
	{
	}
	void init()
	{
		static bool initialized = false;
		if (!initialized) {
			ProjDataInMemory* self = (ProjDataInMemory*)this;
			acqs_templ_.reset(new BinAcquisitionsInMemory
				(self->get_exam_info_sptr(), self->get_proj_data_info_sptr()));
			initialized = true;
		}
	}
	void set_as_template()
	{
		ProjDataInMemory* self = (ProjDataInMemory*)this;
		init();
		acqs_templ_.reset(new BinAcquisitionsInMemory
			(self->get_exam_info_sptr(), self->get_proj_data_info_sptr()));
	}
	PETAcquisitions* same_acquisitions_container()
	{
		ProjDataInMemory* self = (ProjDataInMemory*)this;
		PETAcquisitions* ptr_acq = (PETAcquisitions*)
			new ProjDataInMemory(self->get_exam_info_sptr(),
				self->get_proj_data_info_sptr());
		return ptr_acq;
	}
	boost::shared_ptr<aDataContainer<float> > new_data_container()
	{
		init();
		return boost::shared_ptr<aDataContainer<float> >
			(acqs_templ_->same_acquisitions_container());
	}
	boost::shared_ptr<PETAcquisitions> new_acquisitions_container()
	{
		init();
		return boost::shared_ptr<PETAcquisitions>
			(acqs_templ_->same_acquisitions_container());
	}
};

class BinAcquisitionsInFile : 
	public BinAcquisitions, public ProjDataInterfile {
public:
	BinAcquisitionsInFile(shared_ptr<ExamInfo> const& exam_info_sptr,
		shared_ptr<ProjDataInfo> const& proj_data_info_ptr, 
		const char* filename) :
		ProjDataInterfile(exam_info_sptr, proj_data_info_ptr, filename, 
			std::ios::in | std::ios::out)
	{
	}
	void init()
	{
		static bool initialized = false;
		if (!initialized) {
			ProjDataInMemory* self = (ProjDataInMemory*)this;
			acqs_templ_.reset(new BinAcquisitionsInMemory
			(self->get_exam_info_sptr(), self->get_proj_data_info_sptr()));
			initialized = true;
		}
	}
	void set_as_template()
	{
		ProjDataInterfile* self = (ProjDataInterfile*)this;
		std::string filename = SIRFUtilities::scratch_file_name();
		init();
		acqs_templ_.reset(new BinAcquisitionsInFile
			(self->get_exam_info_sptr(), self->get_proj_data_info_sptr(), 
				filename.c_str()));
	}
	PETAcquisitions* same_acquisitions_container()
	{
		ProjDataInterfile* self = (ProjDataInterfile*)this;
		std::string filename = SIRFUtilities::scratch_file_name();
		PETAcquisitions* ptr_acq = (PETAcquisitions*)
			new ProjDataInterfile(self->get_exam_info_sptr(),
				self->get_proj_data_info_sptr(), filename.c_str(),
				std::ios::in | std::ios::out);
		return ptr_acq;
	}
	boost::shared_ptr<aDataContainer<float> > new_data_container()
	{
		init();
		return boost::shared_ptr<aDataContainer<float> >
			(acqs_templ_->same_acquisitions_container());
	}
	boost::shared_ptr<PETAcquisitions> new_acquisitions_container()
	{
		init();
		return boost::shared_ptr<PETAcquisitions>
			(acqs_templ_->same_acquisitions_container());
	}
};
#endif

class PETImage : public Image3DF, public aDataContainer<float> {
public:
	boost::shared_ptr<aDataContainer<float> > new_data_container()
	{
		Image3DF* self = (Image3DF*)this;
		PETImage* ptr_image = (PETImage*)self->get_empty_copy();
		aDataContainer<float>* ptr_data = (aDataContainer<float>*)ptr_image;
		return boost::shared_ptr<aDataContainer<float> >(ptr_data);
	}
	unsigned int items()
	{
		return 1;
	}
	float norm()
	{
		return 0;
	}
	float dot(PETImage& other)
	{
		return 0;
	}
	void axpby(float a, const PETImage& x, float b, const PETImage& y)
	{

	}
};

template<class Image>
class PETAcquisitionModel {
public:
	void set_projectors(boost::shared_ptr<ProjectorByBinPair> sptr_projectors)
	{
		sptr_projectors_ = sptr_projectors;
	}
	boost::shared_ptr<ProjectorByBinPair> projectors_sptr()
	{
		return sptr_projectors_;
	}
	void set_additive_term(boost::shared_ptr<ProjData> sptr)
	{
		sptr_add_ = sptr;
	}
	void set_background_term(boost::shared_ptr<ProjData> sptr)
	{
		sptr_background_ = sptr;
	}
	boost::shared_ptr<ProjData> additive_term_sptr()
	{
		return sptr_add_;
	}
	void set_normalisation(boost::shared_ptr<BinNormalisation> sptr)
	{
		sptr_normalisation_ = sptr;
	}
	boost::shared_ptr<BinNormalisation> normalisation_sptr()
	{
		return sptr_normalisation_;
	}
	void set_bin_efficiency(boost::shared_ptr<ProjData> sptr_data)
	{
		boost::shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data));
		inv_(sptr.get(), MIN_BIN_EFFICIENCY);
		sptr_normalisation_.reset(new BinNormalisationFromProjData(sptr));
	}
	void set_normalisation(boost::shared_ptr<ProjData> sptr_data)
	{
		boost::shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data));
		sptr_normalisation_.reset(new BinNormalisationFromProjData(sptr));
	}
	void cancel_background_term()
	{
		sptr_background_.reset();
	}
	void cancel_additive_term()
	{
		sptr_add_.reset();
	}
	void cancel_normalisation()
	{
		sptr_normalisation_.reset();
	}

	virtual Succeeded set_up(
		boost::shared_ptr<ProjData> sptr_acq,
		boost::shared_ptr<Image> sptr_image)
	{
		Succeeded s = Succeeded::no;
		if (sptr_projectors_.get()) {
			//writeText("setting up PETAcquisitionModel\n");
			s = sptr_projectors_->set_up
				(sptr_acq->get_proj_data_info_sptr(), sptr_image);
			sptr_acq_template_ = sptr_acq;
			sptr_image_template_ = sptr_image;
		}
		return s;
	}

	boost::shared_ptr<ProjData> forward(const Image& image, const char* file = 0)
	{
		boost::shared_ptr<ProjData> sptr_fd;
		//PETAcquisitions* ptr_templ = (PETAcquisitions*)sptr_acq_template_.get();
		//boost::shared_ptr<PETAcquisitions> sptr_templ = ptr_templ->new_acquisitions_container();
		//sptr_fd = boost::dynamic_pointer_cast<ProjData, PETAcquisitions>(ptr_templ->new_acquisitions_container());

		if (file && strlen(file) > 0)
			sptr_fd.reset(
			new ProjDataInterfile(sptr_acq_template_->get_exam_info_sptr(),
			sptr_acq_template_->get_proj_data_info_sptr(), file, std::ios::in | std::ios::out));
		else
			sptr_fd.reset(
			new ProjDataInMemory(sptr_acq_template_->get_exam_info_sptr(),
			sptr_acq_template_->get_proj_data_info_sptr()));

		sptr_projectors_->get_forward_projector_sptr()->forward_project
			(*sptr_fd, image);
		//if (file && strlen(file) > 0) {
		//	sptr_fd.reset();
		//	sptr_fd = ProjData::read_from_file(file);
		//}

		if (sptr_add_.get()) {
			add_(sptr_fd, sptr_add_);
			std::cout << "additive term added\n";
		}
		else
			std::cout << "no additive term added\n";

		if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
			sptr_normalisation_->undo(*sptr_fd, 0, 1);
			//sptr_normalisation_->apply(*sptr_fd, 0, 1);
			std::cout << "normalisation applied\n";
		}
		else
			std::cout << "no normalisation applied\n";

		if (sptr_background_.get()) {
			add_(sptr_fd, sptr_background_);
			std::cout << "background term added\n";
		}
		else
			std::cout << "no background term added\n";

		return sptr_fd;
	}

	boost::shared_ptr<Image> backward(const ProjData& ad)
	{
		boost::shared_ptr<Image> sptr_im(sptr_image_template_->clone());
		sptr_im->fill(0.0);

		if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
			std::cout << "applying normalisation...\n";
			ProjDataInMemory adc(ad);
			sptr_normalisation_->undo(adc, 0, 1);
			sptr_projectors_->get_back_projector_sptr()->back_project
				(*sptr_im, adc);
		}
		else
			sptr_projectors_->get_back_projector_sptr()->back_project
			(*sptr_im, ad);

		return sptr_im;
	}

protected:
	boost::shared_ptr<ProjectorByBinPair> sptr_projectors_;
	boost::shared_ptr<ProjData> sptr_acq_template_;
	boost::shared_ptr<Image> sptr_image_template_;
	boost::shared_ptr<ProjData> sptr_add_;
	boost::shared_ptr<ProjData> sptr_background_;
	boost::shared_ptr<BinNormalisation> sptr_normalisation_;

private:
	void inv_(ProjData* ptr, float minval)
	{
		size_t size = ptr->size_all();
		float* data = new float[size];
		ptr->copy_to(data);
		inv_(size, data, minval);
		ptr->fill_from(data);
		delete[] data;
	}
	void add_
		(boost::shared_ptr<ProjData> sptr_a, boost::shared_ptr<ProjData> sptr_b)
	{
		size_t size_a = sptr_a->size_all();
		size_t size_b = sptr_b->size_all();
		if (size_a != size_b) {
			std::cout << "ERROR: mismatching sizes " << size_a
				<< " and " << size_b << ", skipping\n";
			return;
		}
		float* data_a = new float[size_a];
		float* data_b = new float[size_b];
		sptr_a->copy_to(data_a);
		sptr_b->copy_to(data_b);
		add_(size_a, data_a, data_b);
		sptr_a->fill_from(data_a);
		delete[] data_a;
		delete[] data_b;
	}
	void inv_(size_t n, float* u, float minval)
	{
		for (size_t i = 0; i < n; i++)
			u[i] = 1 / std::max(minval, u[i]);
	}
	void add_(size_t n, float* u, float* v)
	{
		for (size_t i = 0; i < n; i++)
			u[i] += v[i];
	}
};

template<class Image>
class PETAcquisitionModelUsingMatrix : public PETAcquisitionModel<Image> {
public:
	PETAcquisitionModelUsingMatrix()
	{
		this->sptr_projectors_.reset(new ProjectorPairUsingMatrix);
	}
	void set_matrix(boost::shared_ptr<ProjMatrixByBin> sptr_matrix)
	{
		sptr_matrix_ = sptr_matrix;
		((ProjectorPairUsingMatrix*)this->sptr_projectors_.get())->
			set_proj_matrix_sptr(sptr_matrix);
	}
	boost::shared_ptr<ProjMatrixByBin> matrix_sptr()
	{
		return ((ProjectorPairUsingMatrix*)this->sptr_projectors_.get())->
			get_proj_matrix_sptr();
	}
	virtual Succeeded set_up(
		boost::shared_ptr<ProjData> sptr_acq,
		boost::shared_ptr<Image> sptr_image)
	{
		if (!sptr_matrix_.get())
			return Succeeded::no;
		return PETAcquisitionModel<Image>::set_up(sptr_acq, sptr_image);
	}

private:
	boost::shared_ptr<ProjMatrixByBin> sptr_matrix_;
};

typedef PETAcquisitionModel<Image3DF> AcqMod3DF;
typedef PETAcquisitionModelUsingMatrix<Image3DF> AcqModUsingMatrix3DF;
typedef boost::shared_ptr<AcqMod3DF> sptrAcqMod3DF;

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
	void set_acquisition_model(boost::shared_ptr<AcqMod3DF> sptr)
	{
		sptr_am_ = sptr;
		AcqMod3DF& am = *sptr;
		set_projector_pair_sptr(am.projectors_sptr());
		if (am.additive_term_sptr().get())
			set_additive_proj_data_sptr(am.additive_term_sptr());
		if (am.normalisation_sptr().get())
			set_normalisation_sptr(am.normalisation_sptr());
	}
	boost::shared_ptr<AcqMod3DF> acquisition_model_sptr()
	{
		return sptr_am_;
	}
private:
	boost::shared_ptr<AcqMod3DF> sptr_am_;
};

class xSTIR_IterativeReconstruction3DF :
	public IterativeReconstruction<Image3DF> {
public:
	bool post_process() {
		if (this->output_filename_prefix.length() < 1)
			this->set_output_filename_prefix("reconstructed_image");
		return post_processing();
	}
	Succeeded setup(sptrImage3DF const& image) {
		//std::cout << this->output_filename_prefix << '\n';
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

inline int xSTIR_getImageDimensions(const Image3DF& image, int* dim)
{
	dim[0] = 0;
	dim[1] = 0;
	dim[2] = 0;
	Coordinate3D<int> min_indices;
	Coordinate3D<int> max_indices;
	if (!image.get_regular_range(min_indices, max_indices))
		return -1;
	dim[0] = max_indices[1] - min_indices[1] + 1;
	dim[1] = max_indices[2] - min_indices[2] + 1;
	dim[2] = max_indices[3] - min_indices[3] + 1;
	return 0;
}

#endif