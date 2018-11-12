/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include "dynsim_deformer.h"


class DynSimDeformerTester {

private:
	friend class DynamicSimulationDeformer;

public:
	static bool test_deform_contrast_generator( void );
	static bool test_SIRFImageDataDeformation_memory_behavior( void );
	static bool test_deform_pet_contrast_generator( void ) ;
};

