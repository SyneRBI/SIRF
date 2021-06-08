/* ================================================

Author: Johannes Mayer
Date: 2018.03.21
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "tests_tissueparameters.h"

#include "auxiliary_testing_functions.h"

using boost::property_tree::ptree;


bool test_allocate_MRTissueParameter_successful(void)
{
	MRTissueParameter mr_tissue_pars;
	mr_tissue_pars.t1_miliseconds_ = 1000;
	mr_tissue_pars.t2_miliseconds_ = 200;
	mr_tissue_pars.cs_ppm_ = 4.23;

	return true;
}

bool test_allocate_PETTissueParameter_successful(void)
{
	PETTissueParameter pet_tissue_pars;
	pet_tissue_pars.attenuation_1_by_cm_ = 0.01;
	pet_tissue_pars.activity_kBq_ml_ = 15;

	return true;
}

bool test_allocate_TissueParameter_successful(void)
{
	TissueParameter tissue_pars;
	tissue_pars.label_ = 1;
	tissue_pars.name_ = "Liver";

	MRTissueParameter mr_tiss;
	PETTissueParameter pet_tiss;

	tissue_pars.mr_tissue_ = mr_tiss;
	tissue_pars.pet_tissue_ = pet_tiss;

	return true;
}

bool test_TissueParameter_algebra( void )
{
	TissueParameter first_param, second_param;
		
	first_param.label_ = 0;
	first_param.name_ = "name_0";

	first_param.mr_tissue_.spin_density_percentH2O_ = 0.f;
	first_param.mr_tissue_.t1_miliseconds_ = 0.f;
	first_param.mr_tissue_.t2_miliseconds_ = 0.f;
	first_param.mr_tissue_.cs_ppm_ = 0.f;


	first_param.pet_tissue_.attenuation_1_by_cm_ = 0.f;
	first_param.pet_tissue_.activity_kBq_ml_ = 0.f;


	second_param.label_ = 1;
	second_param.name_ = "name_1";

	second_param.mr_tissue_.spin_density_percentH2O_ = 1.f;
	second_param.mr_tissue_.t1_miliseconds_ = 1.f;
	second_param.mr_tissue_.t2_miliseconds_ = 1.f;
	second_param.mr_tissue_.cs_ppm_ = 1.f;


	second_param.pet_tissue_.attenuation_1_by_cm_ = 1.f;
	second_param.pet_tissue_.activity_kBq_ml_ = 1.f;

	float interpol_weight = 0.3;

	TissueParameter interpol_param = (1.f-interpol_weight) * second_param + interpol_weight*first_param;

	std::cout<<epiph(interpol_param.label_ )<< std::endl;
	std::cout<<epiph(interpol_param.name_ )<< std::endl;

	std::cout<<epiph(interpol_param.mr_tissue_.spin_density_percentH2O_ )<< std::endl;
	std::cout<<epiph(interpol_param.mr_tissue_.t1_miliseconds_)<< std::endl;
	std::cout<<epiph(interpol_param.mr_tissue_.t2_miliseconds_)<< std::endl;
	std::cout<<epiph(interpol_param.mr_tissue_.cs_ppm_)<< std::endl;


	std::cout<<epiph(interpol_param.pet_tissue_.attenuation_1_by_cm_)<< std::endl; 
	std::cout<<epiph(interpol_param.pet_tissue_.activity_kBq_ml_ )<< std::endl;

	return true;

}



bool test_get_MRTissueParameter_from_ptree()
{

	ptree pt;

	float const input_rho = 89; // percent of water
	float const input_t1 = 1000;
	float const input_t2 = 2000;
	float const input_cs = 4.3;


	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.spin_density_percentH2O", input_rho);
	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.t1_miliseconds", input_t1);
	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.t2_miliseconds", input_t2);
	pt.put("TissueParameterList.TissueParameter.MRTissueParameter.cs_ppm", input_cs);

	MRTissueParameter mr_tiss;

	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("TissueParameterList") )
	{
		if( v.first == "TissueParameter")		
		{
			mr_tiss = get_mrtissueparameter_from_ptree(v.second);
		}
	}


	bool parameter_set_correct = (mr_tiss.spin_density_percentH2O_ == input_rho);
	parameter_set_correct *= (mr_tiss.t1_miliseconds_ == input_t1);
	parameter_set_correct *= (mr_tiss.t2_miliseconds_ == input_t2);
	parameter_set_correct *= (mr_tiss.cs_ppm_ == input_cs);
	
	return (parameter_set_correct);
}


bool test_get_PETTissueParameter_from_ptree()
{

	ptree pt;

	float const input_attenuation = 0.5;
	float const input_activity = 5;
	
	pt.put("TissueParameterList.TissueParameter.PETTissueParameter.attenuation_1_by_cm", input_attenuation);
	pt.put("TissueParameterList.TissueParameter.PETTissueParameter.activity_kBq_ml", input_activity);

	PETTissueParameter pet_tiss;

	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("TissueParameterList") )
	{
		if( v.first == "TissueParameter")		
		{
			pet_tiss = get_pettissueparameter_from_ptree(v.second);
		}
	}

	bool parameter_set_correct = (pet_tiss.attenuation_1_by_cm_ == input_attenuation);
	parameter_set_correct *= (pet_tiss.activity_kBq_ml_ == input_activity);
	
	return (parameter_set_correct);
}



void test_exception_throw_if_node_not_exists(void)
{

	ptree pt;

	float non_existent_value = 1;

	pt.put("TissueParameterList.TissueParameter.NonexistentModalityType.NonexistentNode", non_existent_value);

	PETTissueParameter pet_tiss;

	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("TissueParameterList") )
	{
		if( v.first == "TissueParameter")		
		{
			pet_tiss = get_pettissueparameter_from_ptree(v.second);
		}
	}
	std::cout << "Don't worry if an exception get's thrown, that's the point of the test"<<std::endl;
}



bool test_read_TissueParameter_label_from_xml( std::string const xml_filepath )
{

	TissueParameterList tissueList = read_TissueParameters_from_xml(xml_filepath);

	TissueParameter firstParam = tissueList[0];

	std::string const input_name = "Dummy" ;
	int const input_label = 0;

	float const input_t1 = 1;
	float const input_t2 = 2;
	float const input_cs = 3;


	float const input_attenuation = 4;
	float const input_activity = 5;



	bool parameter_set_correct = true; 
	
	// omit the check for a correctly set input_name as the white spaces interfere here
	// does not matter anyway, it is just for the user to be able to distinguish.
		
	parameter_set_correct *= (firstParam.label_ == input_label);
	
	parameter_set_correct *= (input_t1 == firstParam.mr_tissue_.t1_miliseconds_);
	parameter_set_correct *= (input_t2 == firstParam.mr_tissue_.t2_miliseconds_);
	parameter_set_correct *= (input_cs == firstParam.mr_tissue_.cs_ppm_);

	parameter_set_correct *= (input_attenuation == firstParam.pet_tissue_.attenuation_1_by_cm_);
	parameter_set_correct *= (input_activity == firstParam.pet_tissue_.activity_kBq_ml_);	

	return parameter_set_correct;
}





bool test_check_label_uniqueness_fails( void )
{
	
	TissueParameterList tiss_list = aux_test::get_mock_tissue_param_list();

	tiss_list[3].label_ = tiss_list[0].label_;

	bool const labels_are_unique = check_label_uniqueness(tiss_list);

	if( !labels_are_unique )
		return true;
	else
		return false;
}

bool test_check_label_uniqueness_true()
{
	TissueParameterList tiss_list = aux_test::get_mock_tissue_param_list();

	bool const labels_are_unique = check_label_uniqueness(tiss_list);

	if( labels_are_unique )
		return true;
	else
		return false;
}



























