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

/*!
\file
\ingroup Gadgets Chains Library
\brief Specification file for the library of SIRF generators of xml-definitions
			 of Gadgetron gadget chains.

\author Evgueni Ovtchinnikov
\author CCP PETMR
*/

#ifndef PROCESSOR_CHAINS_LIBRARY
#define PROCESSOR_CHAINS_LIBRARY

#include "cgadgetron_shared_ptr.h"
#include "gadgetron_x.h"

/**
\brief Class for an xml-definition generator for oversampling removing chain.
*/
class RemoveOversamplingProcessor : public AcquisitionsProcessor {
public:
	RemoveOversamplingProcessor()
	{
		gadgetron::shared_ptr<aGadget> sptr_g(new RemoveROOversamplingGadget);
		add_gadget("gadget", sptr_g);
	}
	static const char* class_name()
	{
		return "RemoveOversamplingProcessor";
	}
};

/**
\brief Class for an xml-definition generator for simple fully sampled 
			 reconstruction chain.
*/
class SimpleReconstructionProcessor : public ImagesReconstructor {
public:
	SimpleReconstructionProcessor()
	{
		//class_ = "SimpleReconstructionProcessor";
		gadgetron::shared_ptr<aGadget> sptr_g(new SimpleReconGadgetSet);
		add_gadget("gadget", sptr_g);
	}
	static const char* class_name()
	{
		return "SimpleReconstructionProcessor";
	}
};

/**
\brief Class for an xml-definition generator for simple undersampled
GRAPPA reconstruction chain.
*/
class SimpleGRAPPAReconstructionProcessor : public ImagesReconstructor {
public:
	SimpleGRAPPAReconstructionProcessor()
	{
		//class_ = "SimpleGRAPPAReconstructionProcessor";
		gadgetron::shared_ptr<aGadget> sptr_g1(new AcquisitionAccumulateTriggerGadget);
		gadgetron::shared_ptr<aGadget> sptr_g2(new BucketToBufferGadget);
		gadgetron::shared_ptr<aGadget> sptr_g3
			(new GenericReconCartesianReferencePrepGadget);
		gadgetron::shared_ptr<aGadget> sptr_g4(new GenericReconCartesianGrappaGadget);
		gadgetron::shared_ptr<aGadget> sptr_g5
			(new GenericReconFieldOfViewAdjustmentGadget);
		gadgetron::shared_ptr<aGadget> sptr_g6(new GenericReconImageArrayScalingGadget);
		gadgetron::shared_ptr<aGadget> sptr_g7(new ImageArraySplitGadget);
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

/**
\brief Class for an xml-definition generator for real image extracting chain.
*/
class ExtractRealImagesProcessor : public ImagesProcessor {
public:
	ExtractRealImagesProcessor()
	{
		//class_ = "ExtractRealImagesProcessor";
		gadgetron::shared_ptr<aGadget> sptr_g(new ExtractGadget);
		add_gadget("gadget", sptr_g);
	}
	static const char* class_name()
	{
		return "ExtractRealImagesProcessor";
	}
};

#endif