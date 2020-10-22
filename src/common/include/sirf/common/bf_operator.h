#pragma once

#include "sirf/STIR/stir_x.h"

namespace sirf {
	class BFOperator : public anOperator<STIRImageData> {
	public:
		BFOperator(std::shared_ptr<PETAcquisitionModel> sptr_am) : sptr_am_(sptr_am) {}
		void set_subset(int sub_num)
		{
			sub_num_ = sub_num;
		}
		void set_num_subsets(int num_sub)
		{
			num_sub_ = num_sub;
		}
		virtual std::shared_ptr<STIRImageData> apply(const STIRImageData& image_data) const
		{
			std::shared_ptr<PETAcquisitionData> sptr_fwd = 
				sptr_am_->forward(image_data, sub_num_, num_sub_, true);
			std::shared_ptr<STIRImageData> sptr_bwd = 
				sptr_am_->backward(*sptr_fwd, sub_num_, num_sub_);
			return sptr_bwd;
		}
	private:
		std::shared_ptr<PETAcquisitionModel> sptr_am_;
		int sub_num_ = 0;
		int num_sub_ = 1;
	};
}