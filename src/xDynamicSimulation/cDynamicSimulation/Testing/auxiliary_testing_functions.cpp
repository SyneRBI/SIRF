/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "auxiliary_testing_functions.h"


#include <omp.h>
#include <sstream>


MRTissueParameter aux_test::get_mock_MR_tissue_parameter(void)
{
	MRTissueParameter mr_tissue_pars;
	mr_tissue_pars.spin_density_percentH2O_ = 100;
	mr_tissue_pars.t1_miliseconds_ = 1;
	mr_tissue_pars.t2_miliseconds_ = 2;
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

ISMRMRD::AcquisitionHeader aux_test::get_mock_acquisition_header( void )
{

	ISMRMRD::AcquisitionHeader acq_hdr;

}

std::string aux_test::get_serialized_mock_ismrmrd_header( void )
{

	
	ISMRMRD::IsmrmrdHeader hdr = get_mock_ismrmrd_header();

	std::ostringstream out;

	ISMRMRD::serialize(hdr, out);

	return out.str();
}



ISMRMRD::IsmrmrdHeader aux_test::get_mock_ismrmrd_header( void )
{
	using namespace ISMRMRD;

	IsmrmrdHeader hdr;

	SequenceParameters seq_pars = get_mock_sequence_parameters();
	AcquisitionSystemInformation asi = get_mock_acquisition_system_information();

	// necessary
	hdr.experimentalConditions = get_mock_experimental_conditions();
	hdr.encoding = aux_test::get_mock_encoding_vector();

	// optional 
	hdr.sequenceParameters = Optional<SequenceParameters>(seq_pars); 
	hdr.acquisitionSystemInformation = Optional<AcquisitionSystemInformation>(asi);





	return hdr;

}

ISMRMRD::AcquisitionSystemInformation aux_test::get_mock_acquisition_system_information( void )
{
	ISMRMRD::AcquisitionSystemInformation asi;

	float const field_strength_t = 1.00; 


	asi.systemFieldStrength_T = ISMRMRD::Optional<float>(field_strength_t);
	return asi;

}


ISMRMRD::SequenceParameters aux_test::get_mock_sequence_parameters( void )
{
	
	
	using namespace ISMRMRD;

	typedef std::vector<float> ParType;
	
	SequenceParameters seq_pars;

	ParType TR = {2};
	ParType TE = {1};
	ParType TI = {1};
	ParType flipAngle_deg = {90};
	std::string sequ_type = {"Flash"};
	ParType dE = {0};

	seq_pars.TR = Optional< ParType >(TR);
	seq_pars.TE = Optional< ParType >(TE);
	seq_pars.TI = Optional< ParType >(TI);
	seq_pars.flipAngle_deg = Optional< ParType >(flipAngle_deg);
	seq_pars.sequence_type = Optional< std::string >(sequ_type);
	seq_pars.echo_spacing = Optional< ParType >(dE);

	return seq_pars;

}

ISMRMRD::ExperimentalConditions aux_test::get_mock_experimental_conditions( void )
{
	ISMRMRD::ExperimentalConditions e_con;
	e_con.H1resonanceFrequency_Hz = 42580000;
	return e_con;
}

ISMRMRD::EncodingSpace aux_test::get_mock_encoded_space( void )
{

	ISMRMRD::MatrixSize mat_size( MOCK_DATA_RO_OVERSAMPLING * MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE);
	ISMRMRD::FieldOfView_mm fov;

	float resolution_mm = 1.5f;

	fov.x = 1/resolution_mm;
	fov.y = 1/resolution_mm;
	fov.z = 1/resolution_mm;


	ISMRMRD::EncodingSpace enc_spac;
	enc_spac.matrixSize = mat_size;
	enc_spac.fieldOfView_mm = fov;

	return enc_spac;
}

ISMRMRD::EncodingSpace aux_test::get_mock_recon_space( void )
{

	ISMRMRD::MatrixSize mat_size(MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE);
	ISMRMRD::FieldOfView_mm fov;

	float resolution_mm = 1.5;
	fov.x = resolution_mm * MOCK_DATA_MATRIX_SIZE;
	fov.y = resolution_mm * MOCK_DATA_MATRIX_SIZE;
	fov.z = resolution_mm * MOCK_DATA_MATRIX_SIZE;

	ISMRMRD::EncodingSpace enc_spac;
	enc_spac.matrixSize = mat_size;
	enc_spac.fieldOfView_mm = fov;

	return enc_spac;

}






ISMRMRD::EncodingLimits aux_test::get_mock_encoding_limits( void )
{
	unsigned short const max_PE1 = MOCK_DATA_MATRIX_SIZE;
	unsigned short const max_PE2 = MOCK_DATA_MATRIX_SIZE;
	
	unsigned short const center_PE1 = MOCK_DATA_MATRIX_SIZE/2 - 1;
	unsigned short const center_PE2 = MOCK_DATA_MATRIX_SIZE/2 - 1;

	ISMRMRD::Limit limit_PE1(0, max_PE1, center_PE1);
	ISMRMRD::Limit limit_PE2(0, max_PE2, center_PE2);

	ISMRMRD::EncodingLimits enc_lim;
	enc_lim.kspace_encoding_step_1 = ISMRMRD::Optional<ISMRMRD::Limit>( limit_PE1);
	enc_lim.kspace_encoding_step_2 = ISMRMRD::Optional<ISMRMRD::Limit>( limit_PE2);

	return enc_lim;
}



std::vector< ISMRMRD::Encoding > aux_test::get_mock_encoding_vector( void )
{
	ISMRMRD::Encoding enc;

	enc.trajectory = "Cartesian";

	enc.encodedSpace = get_mock_encoded_space();
	enc.reconSpace = get_mock_recon_space();
	enc.encodingLimits = get_mock_encoding_limits();


	std::vector< ISMRMRD::Encoding > enc_vec;
	enc_vec.push_back( enc );
	return enc_vec;
}



TissueParameterList aux_test::get_mock_tissue_param_list( void )
{
	TissueParameter par1, par2, par3, par4;
	par1.name_ = "fake_one";
	par1.label_ = 1;

	par2.name_ = "fake_two";
	par2.label_ = 2;

	par3.name_ = "fake_three";
	par3.label_ = 3;

	par4.name_ = "fake_four";
	par4.label_ = 4;

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
			labels_list(i) = 1;
		else
			labels_list(i) = 2;
	}

	return labels_list;	
}



ISMRMRD::NDArray<complex_float_t> aux_test::get_mock_ndarray_with_cube( void )
{



	size_t const Nx = MOCK_DATA_MATRIX_SIZE;
	size_t const Ny = MOCK_DATA_MATRIX_SIZE;
	size_t const Nz = MOCK_DATA_MATRIX_SIZE;

	std::vector< size_t > mock_dims;
	mock_dims.push_back(Nx);
	mock_dims.push_back(Ny);
	mock_dims.push_back(Nz);

	ISMRMRD::NDArray<complex_float_t> mock_arr;
	mock_arr.resize(mock_dims);


	float const val = 1;
	complex_float_t const cube_value = std::complex<float>( val, val);

	float const cube_radius = Nx/4;
	std::vector<float> cube_center = {Nx/2, Ny/2, Nz/2};
	//#pragma omp parallel
	for( size_t nz=0; nz<Nz; nz++)
	{
		bool z_is_in_cube = (std::abs(nz - cube_center[2]) < cube_radius);
		for( size_t ny=0; ny<Ny; ny++)
		{
			bool y_is_in_cube = (std::abs(ny - cube_center[1]) < cube_radius);

			for( size_t nx=0; nx<Nx; nx++)
			{
				bool x_is_in_cube = (std::abs(nx - cube_center[0]) < cube_radius);
				if( x_is_in_cube && y_is_in_cube && z_is_in_cube )
					mock_arr(nx,ny,nz) = cube_value;

				else
					mock_arr(nx,ny,nz) = 0;
			}
		}
	}

	return mock_arr;
}


ISMRMRD::Image< complex_float_t > aux_test::get_mock_ismrmrd_image_with_cube( void )
{

	ISMRMRD::NDArray<complex_float_t> mock_arr = get_mock_ndarray_with_cube();
	size_t const *  img_dims = mock_arr.getDims();

	ISMRMRD::Image< complex_float_t > mock_img(img_dims[0], img_dims[1], img_dims[2], 1);
		
	size_t num_elements = mock_img.getNumberOfDataElements();

	for( size_t i=0; i<num_elements; i++)
	{
		*(mock_img.begin() + i) = *(mock_arr.begin() + i);
	}

	return mock_img;
}



ISMRMRD::NDArray< complex_float_t > aux_test::get_mock_csm( void )
{
	size_t const Nx = MOCK_DATA_MATRIX_SIZE;
	size_t const Ny = MOCK_DATA_MATRIX_SIZE;
	size_t const Nz = MOCK_DATA_MATRIX_SIZE;

	std::vector <size_t> csm_dims;
	csm_dims.push_back( Nx );
	csm_dims.push_back( Ny );
	csm_dims.push_back( Ny );
	csm_dims.push_back( MOCK_DATA_NUM_CHANNELS );

	ISMRMRD::NDArray< complex_float_t > csm(csm_dims);

	for( size_t i=0; i<csm.getNumberOfElements(); i++)
		*(csm.begin() +i) = 0;
	
	for(size_t nc=0; nc<MOCK_DATA_NUM_CHANNELS; nc++)
	{
		for(size_t nz=0; nz<Nz; nz++)
		{
			for(size_t ny; ny<Ny; ny++)
			{
				for(size_t nx; nx<Nx; nx++)
				{
					if( nc%2 == 0 && nz < MOCK_DATA_MATRIX_SIZE/2)
						csm(nx, ny, nz, nc) = 1;	
				
					else if( nc%2 == 1 && nz > MOCK_DATA_MATRIX_SIZE/2 )
						csm(nx, ny, nz, nc) = 1;	

				}
			}
		}
	}

	return csm;
}

CoilDataAsCFImage aux_test::get_mock_coildata_as_cfimage( void )
{

	ISMRMRD::NDArray<complex_float_t> mock_csm = get_mock_csm();
	
	size_t const * dummy_size = mock_csm.getDims();
	std::vector <size_t> data_size;
	for(int i=0; i<ISMRMRD::ISMRMRD_NDARRAY_MAXDIM; i++)
	{
		if(dummy_size[i] > 0)
			data_size.push_back(dummy_size[i]);
		else 
			data_size.push_back( 1 ); 
	}

	CoilDataAsCFImage csm_as_img( data_size[0], data_size[1], data_size[2], data_size[3] );
	csm_as_img.set_data( mock_csm.begin() );

	return csm_as_img;

}


ISMRMRD::Acquisition aux_test::get_mock_ismrmrd_acquisition ( void )
{
	ISMRMRD::Acquisition acq(MOCK_DATA_MATRIX_SIZE, MOCK_DATA_NUM_CHANNELS, 0);
	return acq;
}
