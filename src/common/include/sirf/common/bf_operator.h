#pragma once

#include "sirf/STIR/stir_x.h"

namespace sirf {
	class BFOperator : public anOperator<STIRImageData> {
	public:
		BFOperator(std::shared_ptr<PETAcquisitionModel> sptr_am) : sptr_am_(sptr_am) {}
		virtual std::shared_ptr<STIRImageData> apply(const STIRImageData& image_data) const
		{
			std::shared_ptr<PETAcquisitionData> sptr_fwd = sptr_am_->forward(image_data);
			std::shared_ptr<STIRImageData> sptr_bwd = sptr_am_->backward(*sptr_fwd);
			return sptr_bwd;
		}
	private:
		std::shared_ptr<PETAcquisitionModel> sptr_am_;
	};
}