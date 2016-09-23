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
	void drop_background_term()
	{
		sptr_background_.reset();
	}
	void drop_additive_term()
	{
		sptr_add_.reset();
	}
	void drop_normalisation()
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
		if (sptr_background_.get()) {
			std::cout << "adding background term...\n";
			size_t size = sptr_background_->size_all();
			size_t foreground_size = sptr_fd->size_all();
			if (size != foreground_size) {
				std::cout << "wrong background term size " << size
					<< ", must be " << foreground_size << ", skipping\n";
				return sptr_fd;
			}
			double* fdata = new double[size];
			double* bdata = new double[size];
			sptr_fd->copy_to(fdata);
			sptr_background_->copy_to(bdata);
			add_(size, fdata, bdata);
			sptr_fd->fill_from(fdata);
			delete[] fdata;
			delete[] bdata;
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

	void add_(size_t n, double* u, double* v)
	{
		for (size_t i = 0; i < n; i++)
			u[i] += v[i];
	}
};

#endif