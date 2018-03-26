/* ================================================

Author: Johannes Mayer
Date: 2018.03.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "phantom_input.h"




using namespace H5;

ISMRMRD::NDArray< unsigned int > read_segmentation_from_h5( std::string const h5_filename_with_suffix)
{
	

	const H5std_string name_segmentation_dataset = "segmentation";
	H5File file( h5_filename_with_suffix, H5F_ACC_RDONLY );
	
	//DataSet segmentation_dataset = file.openDataSet( name_segmentation_dataset );
/*

	H5T_class_t type_class = segmentation_dataset.getTypeClass();

	if( type_class == H5T_INTEGER )
	{
		IntType integer_type = segmentation_dataset.getIntType();
		      
          
     	H5std_string order_string;
        H5T_order_t order = integer_type.getOrder( order_string );
     	std::cout << order_string << std::endl;
         
                   
        size_t size = integer_type.getSize();
        std::cout << "Data size is " << size << std::endl;


	}
	else
	{
		throw std::runtime_error("Please give an integer array as input for the segmentation.");
	}
*/




	ISMRMRD::NDArray< unsigned int > segmentation;
	return segmentation;
}