#ifndef PROCESSOR_CHAINS_LIBRARY
#define PROCESSOR_CHAINS_LIBRARY

#include "xgadgetron.h"

class RemoveOversamplingProcessor : public AcquisitionsProcessor {
public:
	RemoveOversamplingProcessor()
	{
		//class_ = "RemoveOversamplingProcessor";
		boost::shared_ptr<aGadget> sptr_g(new RemoveROOversamplingGadget);
		add_gadget("gadget", sptr_g);
	}
	static const char* class_name()
	{
		return "RemoveOversamplingProcessor";
	}
};

class SimpleReconstructionProcessor : public ImagesReconstructor {
public:
	SimpleReconstructionProcessor()
	{
		//class_ = "SimpleReconstructionProcessor";
		boost::shared_ptr<aGadget> sptr_g(new SimpleReconGadgetSet);
		add_gadget("gadget", sptr_g);
	}
	static const char* class_name()
	{
		return "SimpleReconstructionProcessor";
	}
};

class SimpleGRAPPAReconstructionProcessor : public ImagesReconstructor {
public:
	SimpleGRAPPAReconstructionProcessor()
	{
		//class_ = "SimpleGRAPPAReconstructionProcessor";
		boost::shared_ptr<aGadget> sptr_g1(new AcquisitionAccumulateTriggerGadget);
		boost::shared_ptr<aGadget> sptr_g2(new BucketToBufferGadget);
		boost::shared_ptr<aGadget> sptr_g3
			(new GenericReconCartesianReferencePrepGadget);
		boost::shared_ptr<aGadget> sptr_g4(new GenericReconCartesianGrappaGadget);
		boost::shared_ptr<aGadget> sptr_g5
			(new GenericReconFieldOfViewAdjustmentGadget);
		boost::shared_ptr<aGadget> sptr_g6(new GenericReconImageArrayScalingGadget);
		boost::shared_ptr<aGadget> sptr_g7(new ImageArraySplitGadget);
		add_gadget("gadget1", sptr_g1);
		add_gadget("gadget2", sptr_g2);
		add_gadget("gadget3", sptr_g3);
		add_gadget("gadget4", sptr_g4);
		add_gadget("gadget5", sptr_g5);
		add_gadget("gadget6", sptr_g6);
		add_gadget("gadget7", sptr_g7);
	}
	static const char* class_name()
	{
		return "SimpleGRAPPAReconstructionProcessor";
	}
};

class ExtractRealImagesProcessor : public ImagesProcessor {
public:
	ExtractRealImagesProcessor()
	{
		//class_ = "ExtractRealImagesProcessor";
		boost::shared_ptr<aGadget> sptr_g(new ExtractGadget);
		add_gadget("gadget", sptr_g);
	}
	static const char* class_name()
	{
		return "ExtractRealImagesProcessor";
	}
};

#endif