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
	
	std::cout <<"Trying to read xml file: " << xml_filepath << std::endl;
	
	try
	{
		std::ifstream in(xml_filepath);
		in.exceptions( in.failbit );

		using boost::property_tree::ptree;
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

				tiss_list.push_back(tiss_par);
			}
		}

		return tiss_list;
	}
	catch ( const std::ios_base::failure &e)
	{
		std::cout	<< "Failed to open " << xml_filepath << "\n"
					<< "Caught an ios_base::failure" << "\n"
					<< "Explanatory string: " << e.what() << "\n"
					<< "Error code: " << e.code() << std::endl;
		throw;					
		
	}

}


