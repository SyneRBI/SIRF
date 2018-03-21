/* Implementation of TissueParameter

author		johannes mayer
date		20180321
institute	PTB Berlin

*/

#include <stdio.h>
#include <iostream>
#include "tissueparameters.h"


TissueParameterList read_TissueParameters_from_xml(std::string const xml_filepath)
{	
	std::cout << xml_filepath << std::endl;
	std::ifstream in(xml_filepath);

	if ( in.is_open() )
		std::cout<< "file was opened successfully" <<std::endl;
	else
		std::cout<< "failed to open " << xml_filepath << std::endl;

	using boost::property_tree::ptree;
	ptree pt;

	read_xml( in, pt);	


	TissueParameterList tiss_list;

	BOOST_FOREACH( ptree::value_type const& v, pt.get_child("TissueParameter") )
	{
		if( v.first == "TissueParameter")		
		{
			TissueParameter tiss_par;
			tiss_par.label_ = v.second.get< unsigned >( "label");
			tiss_par.name_  = v.second.get< std::string > ("name");

			tiss_list.push_back(tiss_par);
		}
	}

	return tiss_list;
}


