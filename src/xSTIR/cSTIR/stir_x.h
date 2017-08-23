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

#include <stdlib.h>

#include <chrono>
#include <fstream>

#include "data_handle.h"
#include "stir_types.h"

//#define MIN_BIN_EFFICIENCY 1.0e-20f
#define MIN_BIN_EFFICIENCY 1.0e-6f

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
		sprintf(buff, "tmp_%d_%lld", calls, ms);
		return std::string(buff);
	}
};

template <typename T>
class aDataContainer {
public:
	virtual ~aDataContainer() {}
	//virtual aDataContainer<T>* new_data_container() = 0;
	virtual boost::shared_ptr<aDataContainer<T> > new_data_container() = 0;
	virtual unsigned int items() = 0;
	virtual float norm() = 0;
	virtual T dot(const aDataContainer<T>& dc) = 0;
	virtual void mult(T a, const aDataContainer<T>& x) = 0;
	virtual void axpby(
		T a, const aDataContainer<T>& x,
		T b, const aDataContainer<T>& y) = 0;
};

class ProjDataFile : public ProjDataInterfile {
public:
	ProjDataFile(shared_ptr<ExamInfo> const& exam_info_sptr,
		shared_ptr<ProjDataInfo> const& proj_data_info_ptr,
		const std::string& filename, const std::ios::openmode open_mode) :
		ProjDataInterfile(exam_info_sptr, proj_data_info_ptr, filename, open_mode)
	{}
	shared_ptr<std::iostream> sino_stream_sptr()
	{
		return sino_stream;
	}
	void close_stream()
	{
		((std::fstream*)sino_stream.get())->close();
	}
	void clear_stream()
	{
		((std::fstream*)sino_stream.get())->clear();
	}
};

class ProjDataScratchFile {
public:
	ProjDataScratchFile(const ProjData& pd)
	{
		_init(pd);
	}
	~ProjDataScratchFile()
	{
		_data->close_stream();
		int err;
		err = std::remove((_filename + ".hs").c_str());
		if (err)
			std::cout << "deleting " << _filename << ".hs "
			<< "failed, please delete manually" << std::endl;
		err = std::remove((_filename + ".s").c_str());
		if (err)
			std::cout << "deleting " << _filename << ".s "
			<< "failed, please delete manually" << std::endl;
	}
	void clear_stream()
	{
		_data->clear_stream();
	}
	boost::shared_ptr<ProjData> data() 
	{
		return _data;
	}
private:
	std::string _filename;
	boost::shared_ptr<ProjDataFile> _data;
	void _init(const ProjData& pd)
	{
		_filename = SIRFUtilities::scratch_file_name();
		_data.reset(new ProjDataFile(pd.get_exam_info_sptr(),
			pd.get_proj_data_info_sptr(),
			_filename, std::ios::in | std::ios::out | std::ios::trunc));
	}
};

class PETAcquisitionData : public aDataContainer < float > {
public:
	virtual ~PETAcquisitionData() {}
	static void set_storage_scheme(std::string scheme) 
	{
		_init();
		_storage_scheme = scheme; 
	}
	void read_from_file(const char* filename)
	{
		_data = ProjData::read_from_file(filename);
	}
	PETAcquisitionData* same_acquisition_data(const ProjData& pd)
	{
		PETAcquisitionData* ptr_ad = new PETAcquisitionData;
		_init();
		if (_storage_scheme[0] == 'm') {
			ptr_ad->_data = boost::shared_ptr<ProjData>
				(new ProjDataInMemory(pd.get_exam_info_sptr(),
				pd.get_proj_data_info_sptr()));
			return ptr_ad;
		}
		else {
			ProjDataScratchFile* file = new ProjDataScratchFile(pd);
			ptr_ad->_data = file->data();
			ptr_ad->_file.reset(file);
			return ptr_ad;
		}
	}
	boost::shared_ptr<PETAcquisitionData> new_acquisition_data()
	{
		return boost::shared_ptr<PETAcquisitionData>
			(same_acquisition_data(*_data));
	}
	boost::shared_ptr<aDataContainer<float> > new_data_container()
	{
		return boost::shared_ptr<aDataContainer<float> >
			(same_acquisition_data(*_data));
	}

	// ProjData accessor/mutator
	boost::shared_ptr<ProjData> data() { return _data; }
	void set_data(boost::shared_ptr<ProjData> data) 
	{
		_data = data;
	}

	// data import/export
	void fill(float v) { _data->fill(v); }
	void fill(PETAcquisitionData& data)
	{
		boost::shared_ptr<ProjData> sptr = data.data();
		_data->fill(*sptr);
	}
	void fill_from(const float* data) { _data->fill_from(data); }
	void copy_to(float* data) { _data->copy_to(data); }

	// data container methods
	unsigned int items() { return 1; }
	float norm();
	float dot(const aDataContainer<float>& x) 
	{
		return 0;
	}
	void mult(float a, const aDataContainer<float>& x);
	void inv(float a, const aDataContainer<float>& x);
	void axpby(float a, const aDataContainer<float>& x,
		float b, const aDataContainer<float>& y);

	// ProjData methods
	int get_num_tangential_poss()
	{
		return _data->get_num_tangential_poss();
	}
	int get_num_views()
	{
		return _data->get_num_views();
	}
	int get_num_sinograms()
	{
		return _data->get_num_sinograms();
	}
	int get_max_segment_num() const
	{
		return _data->get_max_segment_num();
	}
	SegmentBySinogram<float> 
		get_segment_by_sinogram(const int segment_num) const
	{
		return _data->get_segment_by_sinogram(segment_num);
	}
	SegmentBySinogram<float>
		get_empty_segment_by_sinogram(const int segment_num) const
	{
		return _data->get_empty_segment_by_sinogram(segment_num);
	}
	virtual Succeeded set_segment(const SegmentBySinogram<float>& s)
	{
		return _data->set_segment(s);
	}
	boost::shared_ptr<ExamInfo> get_exam_info_sptr() const
	{
		return _data->get_exam_info_sptr();
	}
	boost::shared_ptr<ProjDataInfo> get_proj_data_info_sptr() const
	{
		return _data->get_proj_data_info_sptr();
	}

	void clear_stream()
	{
		if (_file.get())
			_file->clear_stream();
	}

	// ProjData casts
	operator ProjData&() { return *_data; }
	operator const ProjData&() const { return *_data; }
	operator ProjData*() { return _data.get(); }
	operator const ProjData*() const { return _data.get(); }
	operator boost::shared_ptr<ProjData>() { return _data; }
	operator const boost::shared_ptr<ProjData>() const { return _data; }

protected:
	static std::string _storage_scheme;
	//static boost::shared_ptr<PETAcquisitionData> _template;
	static void _init()
	{
		static bool initialized = false;
		if (!initialized) {
			//_template.reset();
			_storage_scheme = "file";
			initialized = true;
		}
	}
	boost::shared_ptr<ProjData> _data;
	boost::shared_ptr<ProjDataScratchFile> _file;
};

class PETAcquisitionDataInFile : public PETAcquisitionData {
public:
	PETAcquisitionDataInFile(const char* filename)
	{
		_data = ProjData::read_from_file(filename);
	}
	PETAcquisitionDataInFile(const ProjData& pd)
	{
		_file.reset(new ProjDataScratchFile(pd));
		_data = _file->data();
	}
};

class PETAcquisitionDataInMemory : public PETAcquisitionData {
public:
	PETAcquisitionDataInMemory(const ProjData& pd)
	{
		_data = boost::shared_ptr<ProjData>
			(new ProjDataInMemory(pd.get_exam_info_sptr(),
			pd.get_proj_data_info_sptr()));
	}
};

class PETImageData : public aDataContainer<float> {
public:
	PETImageData(){}
	PETImageData(const Image3DF& image)
	{
		_data.reset(image.clone());
	}
	PETImageData(const Voxels3DF& v)
	{
		_data.reset(v.clone());
	}
	PETImageData(const ProjDataInfo& pdi)
	{
		_data.reset(new Voxels3DF(pdi));
	}
	PETImageData(std::auto_ptr<Image3DF> ptr)
	{
		_data = ptr;
	}
	PETImageData(boost::shared_ptr<Image3DF> ptr)
	{
		_data = ptr;
	}
	PETImageData* same_image_data()
	{
		PETImageData* ptr_image = new PETImageData;
		ptr_image->_data.reset(_data->get_empty_copy());
		return ptr_image;
	}
	boost::shared_ptr<PETImageData> new_image_data()
	{
		return boost::shared_ptr<PETImageData>(same_image_data());
	}
	boost::shared_ptr<aDataContainer<float> > new_data_container()
	{
		return boost::shared_ptr<aDataContainer<float> >(same_image_data());
	}
	unsigned int items()
	{
		return 1;
	}
	float norm();
	float dot(const aDataContainer<float>& other);
	void mult(float a, const aDataContainer<float>& x);
	void axpby(float a, const aDataContainer<float>& x,
		float b, const aDataContainer<float>& y);
	Image3DF& data()
	{
		return *_data;
	}
	const Image3DF& data() const
	{
		return *_data;
	}
	Image3DF* data_ptr()
	{
		return _data.get();
	}
	const Image3DF* data_ptr() const
	{
		return _data.get();
	}
	boost::shared_ptr<Image3DF> data_sptr()
	{
		return _data;
	}
	void fill(float v)
	{
		_data->fill(v);
	}

protected:
	boost::shared_ptr<Image3DF> _data;
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
	//void set_additive_term(boost::shared_ptr<ProjData> sptr)
	void set_additive_term(boost::shared_ptr<PETAcquisitionData> sptr)
	{
		sptr_add_ = sptr;
	}
	boost::shared_ptr<PETAcquisitionData> additive_term_sptr()
	{
		return sptr_add_;
	}
	//void set_background_term(boost::shared_ptr<ProjData> sptr)
	void set_background_term(boost::shared_ptr<PETAcquisitionData> sptr)
	{
		sptr_background_ = sptr;
	}
	boost::shared_ptr<PETAcquisitionData> background_term_sptr()
	{
		return sptr_background_;
	}
	void set_normalisation(boost::shared_ptr<BinNormalisation> sptr)
	{
		sptr_normalisation_ = sptr;
	}
	boost::shared_ptr<BinNormalisation> normalisation_sptr()
	{
		return sptr_normalisation_;
	}
	//void set_bin_efficiency(boost::shared_ptr<ProjData> sptr_data)
	void set_bin_efficiency(boost::shared_ptr<PETAcquisitionData> sptr_data)
	{
		////boost::shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data));
		//boost::shared_ptr<ProjData> sptr(new ProjDataInMemory(*sptr_data->data()));
		//inv_(sptr.get(), MIN_BIN_EFFICIENCY);
		//sptr_normalisation_.reset(new BinNormalisationFromProjData(sptr));

		boost::shared_ptr<PETAcquisitionData>
			sptr_ad(sptr_data->new_acquisition_data());
		sptr_ad->inv(MIN_BIN_EFFICIENCY, *sptr_data);
		sptr_normalisation_.reset
			(new BinNormalisationFromProjData(sptr_ad->data()));
		sptr_normalisation_->set_up(sptr_ad->get_proj_data_info_sptr());

		//std::cout << sptr_normalisation_.get() << std::endl;
		//std::cout << sptr_normalisation_->is_trivial() << std::endl;
		//std::cout << normalisation_sptr()->is_trivial() << std::endl;
		//sptr_ad->clear_stream();
		sptr_norm_ = sptr_ad;
	}
	//void set_normalisation(boost::shared_ptr<ProjData> sptr_data)
	void set_normalisation(boost::shared_ptr<PETAcquisitionData> sptr_data)
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
	void clear_stream()
	{
		if (sptr_norm_.get())
			sptr_norm_->clear_stream();
	}

	virtual Succeeded set_up(
		//boost::shared_ptr<ProjData> sptr_acq,
		boost::shared_ptr<PETAcquisitionData> sptr_acq,
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

	//boost::shared_ptr<ProjData>
	boost::shared_ptr<PETAcquisitionData>
		forward(const Image& image, const char* file = 0)
	{
		boost::shared_ptr<PETAcquisitionData> sptr_ad;
		sptr_ad = sptr_acq_template_->new_acquisition_data();

		boost::shared_ptr<ProjData> sptr_fd = sptr_ad->data();

		sptr_projectors_->get_forward_projector_sptr()->forward_project
			(*sptr_fd, image);

		if (sptr_add_.get()) {
			std::cout << "additive term added...";
			sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_add_);
			//add_(sptr_fd, sptr_add_->data());
			std::cout << "ok\n";
		}
		else
			std::cout << "no additive term added\n";

		clear_stream();
		if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
			std::cout << "normalisation applied...";
			sptr_normalisation_->undo(*sptr_fd, 0, 1);
			//sptr_normalisation_->apply(*sptr_fd, 0, 1);
			std::cout << "ok\n";
		}
		else
			std::cout << "no normalisation applied\n";

		if (sptr_background_.get()) {
			std::cout << "background term added...";
			sptr_ad->axpby(1.0, *sptr_ad, 1.0, *sptr_background_);
			//add_(sptr_fd, sptr_background_->data());
			std::cout << "ok\n";
		}
		else
			std::cout << "no background term added\n";

		//sptr_ad->set_data(sptr_fd);
		//return sptr_fd;
		return sptr_ad;
	}

	boost::shared_ptr<Image> backward(const ProjData& ad)
	{
		boost::shared_ptr<Image> sptr_im(sptr_image_template_->clone());
		sptr_im->fill(0.0);

		if (sptr_normalisation_.get() && !sptr_normalisation_->is_trivial()) {
			std::cout << "applying normalisation...";
			ProjDataInMemory adc(ad);
			std::cout << "ok\n";
			sptr_normalisation_->undo(adc, 0, 1);
			std::cout << "backprojecting...";
			sptr_projectors_->get_back_projector_sptr()->back_project
				(*sptr_im, adc);
			std::cout << "ok\n";
		}
		else
			sptr_projectors_->get_back_projector_sptr()->back_project
			(*sptr_im, ad);

		return sptr_im;
	}

protected:
	boost::shared_ptr<ProjectorByBinPair> sptr_projectors_;
	boost::shared_ptr<PETAcquisitionData> sptr_acq_template_;
	//boost::shared_ptr<ProjData> sptr_acq_template_;
	boost::shared_ptr<Image> sptr_image_template_;
	boost::shared_ptr<PETAcquisitionData> sptr_add_;
	boost::shared_ptr<PETAcquisitionData> sptr_background_;
	//boost::shared_ptr<ProjData> sptr_add_;
	//boost::shared_ptr<ProjData> sptr_background_;
	boost::shared_ptr<BinNormalisation> sptr_normalisation_;
	boost::shared_ptr<PETAcquisitionData> sptr_norm_;
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
		//boost::shared_ptr<ProjData> sptr_acq,
		boost::shared_ptr<PETAcquisitionData> sptr_acq,
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
			set_additive_proj_data_sptr(am.additive_term_sptr()->data());
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