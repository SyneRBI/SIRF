/* ================================================

Author: Johannes Mayer
Date: 2018.04.05
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "sirf/cDynamicSimulation/auxiliary_input_output.h"


#include <ismrmrd/dataset.h>
#include <ismrmrd/ismrmrd.h>

#include "sirf/cSTIR/stir_types.h"

using namespace sirf;


// ++++++++++++++++++ print IO +++++++++++++
void print_io::print_voxelized_geometrical_info( const sirf::ImageData& im )
{

    std::shared_ptr<const VoxelisedGeometricalInfo3D > sptr_geo_info = im.get_geom_info_sptr();

    VoxelisedGeometricalInfo3D::Offset 				curr_offset 	=	sptr_geo_info->get_offset();
    VoxelisedGeometricalInfo3D::Spacing 			curr_spacing 	=	sptr_geo_info->get_spacing();
    VoxelisedGeometricalInfo3D::Size 				curr_size 		=	sptr_geo_info->get_size();
    VoxelisedGeometricalInfo3D::DirectionMatrix 	curr_dirmat 	=	sptr_geo_info->get_direction();


    std::cout << "Offset: (";
    for( int i=0;i<3; i++)
    	std::cout <<  curr_offset[i]<< " / ";
    std::cout << ")"<< std::endl;

    std::cout << "Spacing: (";
    for( int i=0;i<3; i++)
    	std::cout <<  curr_spacing[i]<< " / ";
    std::cout << ")"<< std::endl;

    std::cout << "Size: (";
    for( int i=0;i<3; i++)
    	std::cout <<  curr_size[i]<< " / ";
    std::cout << ")"<< std::endl;

    std::cout << "Dir mat: \n";
	for( int i=0;i<3; i++)
	for( int j=0;j<3; j++)
	{
		std::cout << curr_dirmat[i][j] << "  ";
		if(j==2)
			std::cout << "\n";
	}
	std::cout << "\n";

}

// ++++++++++++++++++ data_IO ++++++++++++++++

SignalContainer data_io::read_surrogate_signal( const std::string& filename_time_axis, const std::string& filename_signal )
{
	std::vector< TimeAxisType > time_points = read_single_column_txt<TimeAxisType>(filename_time_axis);
	std::vector< SignalAxisType > signal_points = read_single_column_txt<SignalAxisType>(filename_signal);
  	
	SignalContainer signal;

  	if( time_points.size() == signal_points.size())
  	{
  		std::cout << time_points.size() << " signal points read from file." <<std::endl;
  		for( size_t i=0; i<time_points.size(); i++)
  		{
	  		SignalPoint sp (time_points[i], signal_points[i]);
	  		signal.push_back( sp );
	  	}
  	}
  	else
  		throw std::runtime_error( "The two files given dont have the same number of data points in them." );

  	return signal;
}

void data_io::write_PET_image_to_hv( const std::string& filename_without_ext,const sirf::STIRImageData& img)
{
	std::stringstream stream_filename; 
	stream_filename << filename_without_ext << ".hv";
	
	std::cout << "Writing PET image ... ";
    img.write(stream_filename.str());
	std::cout << "... finished." << std::endl;
}



// ++++++++++++++++++ mr_IO ++++++++++++++++++

ISMRMRD::IsmrmrdHeader mr_io::read_ismrmrd_header( std::string path_ismrmrd_h5_file_with_ext)
{

	ISMRMRD::Dataset d(path_ismrmrd_h5_file_with_ext.c_str(),"dataset", false);

	std::string xml;
	d.readHeader(xml);

	ISMRMRD::IsmrmrdHeader hdr;
	ISMRMRD::deserialize(xml.c_str(), hdr);
	
	return hdr;

}

AcquisitionsVector mr_io::read_ismrmrd_acquisitions( std::string path_ismrmrd_h5_file_with_ext )
{

	ISMRMRD::IsmrmrdHeader hdr = mr_io::read_ismrmrd_header( path_ismrmrd_h5_file_with_ext );
	std::ostringstream out;
	serialize(hdr, out);

	std::cout<< "Started reading acquisitions from " << path_ismrmrd_h5_file_with_ext << std::endl;

	ISMRMRD::Dataset d(path_ismrmrd_h5_file_with_ext.c_str(),"dataset", false);

	uint32_t num_acquis = d.getNumberOfAcquisitions();

	AcquisitionsVector acq_vec(out.str());

	for( uint32_t i_acqu=0; i_acqu<num_acquis; i_acqu++)
	{
		ISMRMRD::Acquisition acq;

		if ((i_acqu%1000) == 0 )
			std::cout << float(i_acqu)/num_acquis*100.f << " % " << std::endl;

		d.readAcquisition( i_acqu, acq);


		if( acq.isFlagSet( ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT ))
		{
			std::cout << "Acquisition # " << i_acqu <<" omitted due to it being a noise sample." <<std::endl;
			continue;
		}

		acq_vec.append_acquisition( acq );

	}

	std::cout<< "Finished reading acquisitions from " << path_ismrmrd_h5_file_with_ext << std::endl;

	return acq_vec;

}



// ++++++++++++++++++ pet_IO ++++++++++++++++++