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

class SimpleGRAPPAReconstructionProcessor : public ImagesReconstructor {
public:
	SimpleGRAPPAReconstructionProcessor()
	{
		boost::shared_ptr<aGadget> sptr_g1(new AcqAccTrigGadget);
		boost::shared_ptr<aGadget> sptr_g2(new BucketToBuffGadget);
		boost::shared_ptr<aGadget> sptr_g3(new PrepRefGadget);
		boost::shared_ptr<aGadget> sptr_g4(new CartesianGrappaGadget);
		boost::shared_ptr<aGadget> sptr_g5(new FOVAdjustmentGadget);
		boost::shared_ptr<aGadget> sptr_g6(new ScalingGadget);
		boost::shared_ptr<aGadget> sptr_g7(new ImgArrSplitGadget);
		add_gadget("gadget1", sptr_g1);
		add_gadget("gadget2", sptr_g2);
		add_gadget("gadget3", sptr_g3);
		add_gadget("gadget4", sptr_g4);
		add_gadget("gadget5", sptr_g5);
		add_gadget("gadget6", sptr_g6);
		add_gadget("gadget7", sptr_g7);
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