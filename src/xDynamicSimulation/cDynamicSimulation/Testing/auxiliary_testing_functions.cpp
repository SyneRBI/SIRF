/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "auxiliary_testing_functions.h"




MRTissueParameter aux_test::get_mock_MR_tissue_parameter(void)
{
	MRTissueParameter mr_tissue_pars;
	mr_tissue_pars.spin_density_percentH2O_ = 100;
	mr_tissue_pars.t1_miliseconds_ = 0;
	mr_tissue_pars.t2_miliseconds_ = 1;
	mr_tissue_pars.cs_ppm_ = 1;

	return mr_tissue_pars;
}

PETTissueParameter aux_test::get_mock_PET_tissue_parameter(void)
{
	PETTissueParameter pet_tissue_pars;
	pet_tissue_pars.attenuation_1_by_mm_ = 0.01;
	pet_tissue_pars.suv_ = 15;


	return pet_tissue_pars;
}



TissueParameter aux_test::get_mock_tissue_parameter( void )
{

	TissueParameter tiss_par;
	tiss_par.name_ = "mocktissue";
	tiss_par.label_ = 0;

	tiss_par.mr_tissue_ = get_mock_MR_tissue_parameter();
	tiss_par.pet_tissue_ = get_mock_PET_tissue_parameter();
	return tiss_par;
}


TissueParameterList aux_test::get_mock_tissue_param_list( void )
{
	TissueParameter par1, par2, par3, par4;
	par1.name_ = "fake_one";
	par1.label_ = 0;


	par2.name_ = "fake_two";
	par2.label_ = 1;

	par3.name_ = "fake_three";
	par3.label_ = 2;

	par4.name_ = "fake_four";
	par4.label_ = 3;

	TissueParameterList tiss_list;
	
	tiss_list.push_back(par1);
	tiss_list.push_back(par2);
	tiss_list.push_back(par3);
	tiss_list.push_back(par4);

	return tiss_list;	
}

LabelArray aux_test::get_mock_label_array( void )
{
	
	std::vector< size_t > labels_dims = {2,2,2};
	LabelArray labels_list(labels_dims);

	for( int i=0; i< labels_list.getNumberOfElements(); i++)
	{
		if( i< labels_list.getNumberOfElements()/2 )
			labels_list(i) = 0;
		else
			labels_list(i) = 1;
	}

	return labels_list;	
}

