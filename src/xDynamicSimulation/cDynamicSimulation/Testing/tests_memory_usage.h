/* ================================================

Author: Johannes Mayer
Date: 2018.11.26
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#pragma once



#include <memory>
#include <vector>
#include <stdexcept>

#include "gadgetron_data_containers.h"
#include <ismrmrd/ismrmrd.h>


class ObjectAcquisVect{

public:
	ObjectAcquisVect(){};

	ISMRMRD::Acquisition& get_acq( size_t num )
	{
		if( num < vec_acqs_.size())
			return vec_acqs_[num];
		else
		{	
			throw std::runtime_error("Try to access larger index than present.");
		}
	};
	
	std::shared_ptr< ISMRMRD::Acquisition > get_acq_sptr(size_t num)
	{
		if( num < vec_acqs_.size())
			return std::shared_ptr< ISMRMRD::Acquisition > (&vec_acqs_[num]) ;
		else
		{	
			throw std::runtime_error("Try to access larger index than present.");
		}

	};

	void append_acq( const ISMRMRD::Acquisition& acq )
	{
		vec_acqs_.push_back( acq );
	};

protected:
	std::vector< ISMRMRD::Acquisition > vec_acqs_;

};

class PtrAcquisVect{

public:
	std::shared_ptr< ISMRMRD::Acquisition > get_sptr_acq( size_t num )
	{
		if( num < vec_sptr_acqs_.size() ) 
			return vec_sptr_acqs_[num];
		else
		{
			throw std::runtime_error("Try to access larger index than present.");
		}
		
	};

	void append_sptr_acq( std::shared_ptr< ISMRMRD::Acquisition > sptr_acq )
	{
		vec_sptr_acqs_.push_back( sptr_acq );
	};

protected:
	std::vector< std::shared_ptr< ISMRMRD::Acquisition > > vec_sptr_acqs_;

};


namespace tests_memory{


bool test_acquisition_memory( void );
bool test_acquisition_vector_memory( void );

}