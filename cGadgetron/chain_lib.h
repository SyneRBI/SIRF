#ifndef PROCESSOR_CHAINS_LIBRARY
#define PROCESSOR_CHAINS_LIBRARY

#include "xgadgetron.h"

class RemoveOversamplingProcessor : public AcquisitionsProcessor {
public:
	RemoveOversamplingProcessor()
	{
		boost::shared_ptr<aGadget> sptr_g(new RemoveOversamplingGadget);
		add_gadget("gadget", sptr_g);
	}
};

class SimpleReconstructionProcessor : public ImagesReconstructor {
public:
	SimpleReconstructionProcessor()
	{
		boost::shared_ptr<aGadget> sptr_g(new SimpleReconstructionGadgetSet);
		add_gadget("gadget", sptr_g);
	}
};

class ExtractRealImagesProcessor : public ImagesProcessor {
public:
	ExtractRealImagesProcessor()
	{
		boost::shared_ptr<aGadget> sptr_g(new ExtGadget);
		add_gadget("gadget", sptr_g);
	}
};

#endif