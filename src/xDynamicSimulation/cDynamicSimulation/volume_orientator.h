/* ================================================

Author: Johannes Mayer
Date: 2018.09.18
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#pragma once

#include "gadgetron_data_containers.h"
#include "localised_exception.h"



namespace sirf{

enum ReadoutDir {ro_dir_x, ro_dir_y, ro_dir_z};

class aVolumeOrientator{


public:

	aVolumeOrientator(){};
	aVolumeOrientator(ReadoutDir dir):readout_direction_wrt_input_(dir){};

	void set_readout_direction( ReadoutDir direction){ this->readout_direction_wrt_input_ = direction;};
	ReadoutDir get_readout_direction( ReadoutDir direction){ return this->readout_direction_wrt_input_;};

	template <typename T>
	ISMRMRD::Image<T> reorient_image( ISMRMRD::Image<T> img)
	{

	    uint16_t Nx = img.getMatrixSizeX();
	    uint16_t Ny = img.getMatrixSizeY();
	    uint16_t Nz = img.getMatrixSizeZ(); 
	    uint16_t Nc = img.getNumberOfChannels();

	    std::vector< uint16_t > vol_sizes{Nx, Ny, Nz};
	    vol_sizes = this->reshuffle_indices<uint16_t>( vol_sizes );

	    ISMRMRD::Image<T> reoriented_volume( img );
	    reoriented_volume.resize(vol_sizes[0], vol_sizes[1], vol_sizes[2], Nc);

	    std::vector<int> lin_index{0,1,2};
	    auto shuffled_index = this->reshuffle_indices( lin_index );

	    std::vector< size_t > offsets_for_inversion{0,0,0};
	    std::vector< int > sign_for_inversion{1,1,1};

	    // offsets_for_inversion[ (readout_direction_wrt_input_ - 1 + 3)%3  ] = vol_sizes[(readout_direction_wrt_input_ - 1 + 3)%3]-1; 
	    // offsets_for_inversion[ (readout_direction_wrt_input_ - 2 + 3)%3  ] = vol_sizes[(readout_direction_wrt_input_ - 1 + 3)%3]-1; 
	    // sign_for_inversion[ (readout_direction_wrt_input_ - 1 + 3)%3  ] = -1;
	    // sign_for_inversion[ (readout_direction_wrt_input_ - 2 + 3)%3  ] = -1;

	    for( size_t nc=0; nc<Nc; nc++)
	    for( size_t nz=0; nz<Nz; nz++)
	    for( size_t ny=0; ny<Ny; ny++)
	    for( size_t nx=0; nx<Nx; nx++){

	    	std::vector< size_t > vol_access{nx,ny,nz};

	    	size_t const access_x = offsets_for_inversion[0] + sign_for_inversion[0]*vol_access[ shuffled_index[0] ];
	    	size_t const access_y = offsets_for_inversion[1] + sign_for_inversion[1]*vol_access[ shuffled_index[1] ];
	    	size_t const access_z = offsets_for_inversion[2] + sign_for_inversion[2]*vol_access[ shuffled_index[2] ];

	    	reoriented_volume( access_x, access_y, access_z ,nc) = 	img(nx, ny, nz, nc);
	    }

	    return reoriented_volume;
	}
	

protected:
	
	ReadoutDir readout_direction_wrt_input_;
	
	int sign_permutation;

	template< typename T >
	std::vector<T> reshuffle_indices( std::vector<T> &indices)
	{
		if(indices.size()> 3)
			throw LocalisedException("Only give 3 indices to shuffle volume orientation.", __FILE__, __LINE__);

		if(readout_direction_wrt_input_ > 2)
			throw LocalisedException("Please give only readout directions 0 for x, 1 for y or 2 for z..", __FILE__, __LINE__);

		std::vector<T> shuffled_idxs;
		for(int i=0; i<3; i++)
		{
			shuffled_idxs.push_back( indices[ (i +  this->readout_direction_wrt_input_)%3]);
		}
		return shuffled_idxs;
	}
};

}