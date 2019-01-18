/* Implementation of TissueParameter

author		johannes mayer
date		20180321
institute	PTB Berlin

*/

#include <stdio.h>
#include <iostream>
#include "sirf/cDynamicSimulation/tissueparameters.h"

using boost::property_tree::ptree;



MRTissueParameter operator* (float  const x, const MRTissueParameter& a_mr) 
{
	MRTissueParameter mr_tiss;
	mr_tiss.spin_density_percentH2O_ = x * a_mr.spin_density_percentH2O_;
	mr_tiss.t1_miliseconds_ = x * a_mr.t1_miliseconds_;
	mr_tiss.t2_miliseconds_ = x * a_mr.t2_miliseconds_;
	mr_tiss.cs_ppm_ = x * a_mr.cs_ppm_;
	
	return mr_tiss;
}

PETTissueParameter operator* (float const x, const PETTissueParameter& a_pet) 
{
	PETTissueParameter pet_tiss;
	pet_tiss.attenuation_1_by_cm_ = x * a_pet.attenuation_1_by_cm_;
	pet_tiss.activity_kBq_ml_ = x * a_pet.activity_kBq_ml_;
	
	return pet_tiss;
}

TissueParameter operator* (float const  x, const TissueParameter& a_tiss) 
{ 
	TissueParameter tiss;
	tiss.label_ = a_tiss.label_;
	tiss.name_ =  a_tiss.name_;
	
	tiss.mr_tissue_ =  x * a_tiss.mr_tissue_;
	tiss.pet_tissue_ =  x * a_tiss.pet_tissue_;

	return tiss;
}





TissueParameterList read_TissueParameters_from_xml(const std::string& xml_filepath)
{	
	
	std::cout <<"Reading xml file: " << xml_filepath << "." << std::endl;
	
	try
	{
		std::ifstream in(xml_filepath);
		in.exceptions( in.failbit );


		ptree pt;

		read_xml( in, pt);	

		TissueParameterList tiss_list;

		BOOST_FOREACH( ptree::value_type const& v, pt.get_child("TissueParameterList") )
		{

			
			if( v.first == "TissueParameter")		
			{	
				TissueParameter tiss_par;
				tiss_par.label_ = v.second.get< unsigned >( "label");
				tiss_par.name_  = v.second.get< std::string > ("name");

				
				tiss_par.mr_tissue_ = get_mrtissueparameter_from_ptree( v.second );
				tiss_par.pet_tissue_ = get_pettissueparameter_from_ptree( v.second );
				tiss_list.push_back(tiss_par);
			}
		}

		if(tiss_list.size() == 0)
		{
			throw std::range_error( "The TissueParameterList of your xml file does not contain a single TissueParameter node." );
		}
		if( !check_label_uniqueness(tiss_list) )
		{
			throw std::runtime_error( "The TissueParameterList of your xml file contains multiple parameters with the same label. Please only assign unique labels." );	
		}


		std::cout <<"Reading of: " << xml_filepath << " finished without exception." << std::endl;

		return tiss_list;
	}

	// error handling 

	catch ( const std::ios_base::failure &e)
	{

		std::cout	<< "Failed to open " << xml_filepath << "\n"
		<< "Caught an ios_base::failure" << "\n"
		<< "Explanatory string: " << e.what() << "\n"
		<< "Error code: " << e.code() << std::endl;
		throw;					
	}
	catch ( const boost::property_tree::ptree_bad_path &e)
	{
		std::cout << "Caught bad_path exception " <<std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "You probably forgot to name an essential key for tissue parameters in your xml." << std::endl;

	}
	catch ( const std::range_error &e)
	{
		std::cout << "Caught range error exception " <<std::endl;
		std::cout << e.what() << std::endl;		
	}
}


MRTissueParameter get_mrtissueparameter_from_ptree(boost::property_tree::ptree pt)
{

	MRTissueParameter mr_tiss;

	try
	{
		ptree mr_tissue_tree = pt.get_child("MRTissueParameter");


		mr_tiss.spin_density_percentH2O_ = mr_tissue_tree.get <float> ("spin_density_percentH2O");
		mr_tiss.t1_miliseconds_ = mr_tissue_tree.get <float> ("t1_miliseconds");
		mr_tiss.t2_miliseconds_ = mr_tissue_tree.get <float> ("t2_miliseconds");
		mr_tiss.cs_ppm_ = mr_tissue_tree.get <float> ("cs_ppm");
	}
	catch( const boost::property_tree::ptree_bad_path &e) 
	{	
		std::cout << "Caught bad_path exception " <<std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "You probably forgot to name an essential key for MR contrast in your xml." << std::endl;
	}

	return mr_tiss;
}

PETTissueParameter get_pettissueparameter_from_ptree(boost::property_tree::ptree const pt)
{

	PETTissueParameter pet_tiss;

	try
	{
		ptree pet_tissue_tree = pt.get_child("PETTissueParameter");
		pet_tiss.attenuation_1_by_cm_ = pet_tissue_tree.get <float> ("attenuation_1_by_cm");
		pet_tiss.activity_kBq_ml_ = pet_tissue_tree.get <float> ("activity_kBq_ml");
	}
	catch( const boost::property_tree::ptree_bad_path &e) 
	{	
		std::cout << "Caught bad_path exception." <<std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "You probably forgot to list an essential key for PET contrast in your xml." << std::endl;
	}
	return pet_tiss;	
}

bool check_label_uniqueness( const TissueParameterList&  tiss_list)
{

	size_t const num_tissue_params = tiss_list.size();

	std::vector <LabelType> all_labels;
	all_labels.resize(num_tissue_params);
	
	for (size_t i=0; i<num_tissue_params; i++)
		all_labels[i] = tiss_list[i].label_;

	std::vector<LabelType>::iterator it;

	for (size_t i=0; i<num_tissue_params-1; i++)
	{
		LabelType val_to_find = all_labels[i];
		auto start_iterator = std::next(all_labels.begin(), i+1);

		it = find(start_iterator, all_labels.end(), val_to_find);		

		if( it != all_labels.end() ) 		
			return false;
		
	}

	return true;		
}