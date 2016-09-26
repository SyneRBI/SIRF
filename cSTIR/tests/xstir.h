#ifndef MORE_XSTIR_TYPES
#define MORE_XSTIR_TYPES

template<class Image>
class PETAcquisitionModel {
public:
	PETAcquisitionModel(
		boost::shared_ptr<ProjectorByBinPair> sptr_projectors,
		boost::shared_ptr<ProjData> sptr_acq,
		boost::shared_ptr<Image> sptr_image)
	{
		Succeeded s =
			sptr_projectors->set_up(sptr_acq->get_proj_data_info_sptr(), sptr_image);
		sptr_projectors_ = sptr_projectors;
		sptr_acq_template_ = sptr_acq;
		sptr_image_template_ = sptr_image;
	}
	void set_background_term(boost::shared_ptr<ProjData> sptr)
	{
		sptr_background_ = sptr;
	}
	void set_additive_term(boost::shared_ptr<ProjData> sptr)
	{
		sptr_add_ = sptr;
	}
	void set_normalisation(boost::shared_ptr<BinNormalisation> sptr)
	{
		sptr_normalisation_ = sptr;
	}
	void set_normalisation(boost::shared_ptr<ProjData> sptr)
	{
		sptr_normalisation_.reset(new BinNormalisationFromProjData(sptr));
	}
	//void set_multiplicative_term(boost::shared_ptr<ProjData> sptr)
	//{
	//	sptr_mult_ = sptr;
	//}
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
	boost::shared_ptr<ProjData> forward(const Image& image)
	{
		boost::shared_ptr<ProjData> sptr_fd(
			new ProjDataInMemory(sptr_acq_template_->get_exam_info_sptr(),
			sptr_acq_template_->get_proj_data_info_sptr()));
		sptr_projectors_->get_forward_projector_sptr()->forward_project
			(*sptr_fd, image);
		if (sptr_add_.get()) {
			add_(sptr_fd, sptr_add_);
			std::cout << "adding additive term...\n";
		}
		else
			std::cout << "no additive term to add\n";
		if (sptr_normalisation_.get()) {
			std::cout << "applying normalisation...\n";
			sptr_normalisation_->undo(*sptr_fd, 0, 1);
			//sptr_normalisation_->apply(*sptr_fd, 0, 1);
		}
		if (sptr_background_.get()) {
			add_(sptr_fd, sptr_background_);
			std::cout << "adding background term...\n";
		}
		else
			std::cout << "no background term to add\n";
		return sptr_fd;
	}

private:
	boost::shared_ptr<ProjectorByBinPair> sptr_projectors_;
	boost::shared_ptr<ProjData> sptr_acq_template_;
	boost::shared_ptr<Image> sptr_image_template_;
	//boost::shared_ptr<ProjData> sptr_mult_;
	boost::shared_ptr<ProjData> sptr_add_;
	boost::shared_ptr<ProjData> sptr_background_;
	boost::shared_ptr<BinNormalisation> sptr_normalisation_;

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
		double* data_a = new double[size_a];
		double* data_b = new double[size_b];
		sptr_a->copy_to(data_a);
		sptr_b->copy_to(data_b);
		add_(size_a, data_a, data_b);
		sptr_a->fill_from(data_a);
		delete[] data_a;
		delete[] data_b;
	}
	void add_(size_t n, double* u, double* v)
	{
		for (size_t i = 0; i < n; i++)
			u[i] += v[i];
	}
};

#endif