/* ================================================

Author: Johannes Mayer
Date: 2018.04.06
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once

#include "sirf/cGadgetron/encoding.h"


namespace test_enc
{
	bool test_cube_input();
}

class CartesianEncodingTester
{
public:
	static bool test_sample_fourier_space( void );

};



class RPETester
{

public:
	static bool test_sample_fourier_space( void );

};


class RPETrajectoryPreparationTester{

public:
	
	static bool test_get_set_trajectory();	
	static bool test_get_result_container();	
	
};


class RPESuperInterleavedGoldenCutTester
{
public:
	static bool test_compute_trajectory();

};