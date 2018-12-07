/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "auxiliary_testing_functions.h"


#include <omp.h>
#include <sstream>
#include <math.h>
#include <cmath>
#include <stdexcept>

#include <algorithm>
#include <vector>


using namespace sirf;

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

TissueParameter aux_test::get_mock_tissue_parameter( void )
{

	TissueParameter tiss_par;
	tiss_par.name_ = "mocktissue";
	tiss_par.label_ = 0;

	tiss_par.mr_tissue_ = get_mock_MR_tissue_parameter();
	tiss_par.pet_tissue_ = get_mock_PET_tissue_parameter();
	return tiss_par;
}

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
	pet_tissue_pars.attenuation_1_by_cm_ = 0.01;
	pet_tissue_pars.activity_kBq_ml_ = 15;


	return pet_tissue_pars;
}


std::pair< TissueParameter, TissueParameter> aux_test::get_mock_contrast_signal_extremes( void )
{
	std::pair< TissueParameter, TissueParameter> output;

	TissueParameter	tiss_at_0, tiss_at_1;

	tiss_at_0.name_ = "dynamic_par";
	tiss_at_0.label_ = 10000;

	tiss_at_0.mr_tissue_.spin_density_percentH2O_ = 80;
	tiss_at_0.mr_tissue_.t1_miliseconds_ = 1157;
	tiss_at_0.mr_tissue_.t2_miliseconds_ = 44;
	tiss_at_0.mr_tissue_.cs_ppm_ = 0;
	
	tiss_at_0.pet_tissue_.attenuation_1_by_cm_= 0;
	tiss_at_0.pet_tissue_.activity_kBq_ml_= 0;

	output.first = tiss_at_0;

	tiss_at_1.name_ = "dynamic_par";
	tiss_at_1.label_ = 10000;

	tiss_at_1.mr_tissue_.spin_density_percentH2O_ = 80;
	tiss_at_1.mr_tissue_.t1_miliseconds_ = 600;
	tiss_at_1.mr_tissue_.t2_miliseconds_ = 44;
	tiss_at_1.mr_tissue_.cs_ppm_ = 0;
	
	tiss_at_1.pet_tissue_.attenuation_1_by_cm_= 0;
	tiss_at_1.pet_tissue_.activity_kBq_ml_= 0;

	output.second = tiss_at_1;

	return output;
}

MRContrastGenerator aux_test::get_mock_mr_contrast_generator( void )
{
	LabelArray label_list = get_mock_label_array();

	MRContrastGenerator mr_cont(label_list, XML_TEST_PATH);

	ISMRMRD::IsmrmrdHeader hdr = get_mock_ismrmrd_header();

	mr_cont.set_rawdata_header( hdr);

	return mr_cont;

}

PETContrastGenerator aux_test::get_mock_pet_contrast_generator( void )
{
	LabelArray segmentation_labels = read_segmentation_from_h5( H5_XCAT_PHANTOM_PATH );
	PETContrastGenerator pet_cont_gen( segmentation_labels, XML_XCAT_PATH);

	pet_cont_gen.set_template_image_from_file( PET_TEMPLATE_CONTRAST_IMAGE_DATA_PATH );

	return pet_cont_gen;
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

std::string aux_test::get_serialized_mock_ismrmrd_header( void )
{

	
	ISMRMRD::IsmrmrdHeader hdr = get_mock_ismrmrd_header();

	std::ostringstream out;

	ISMRMRD::serialize(hdr, out);

	return out.str();
}


ISMRMRD::AcquisitionSystemInformation aux_test::get_mock_acquisition_system_information( void )
{
	ISMRMRD::AcquisitionSystemInformation asi;

	asi.receiverChannels = ISMRMRD::Optional<unsigned short > ( MOCK_DATA_NUM_CHANNELS );
	asi.systemFieldStrength_T = ISMRMRD::Optional<float>( MOCK_FIELD_STRENGTH );
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
	ParType dE = {2};

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

ISMRMRD::EncodingSpace aux_test::get_mock_encoded_space( void )
{

	ISMRMRD::MatrixSize mat_size( MOCK_DATA_RO_OVERSAMPLING * MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE,MOCK_DATA_MATRIX_SIZE);
	ISMRMRD::FieldOfView_mm fov;

	float const resolution_mm = MOCK_FOV/MOCK_DATA_MATRIX_SIZE;

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

	fov.x = MOCK_FOV;
	fov.y = MOCK_FOV;
	fov.z = MOCK_FOV;

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


	mock_img.setMatrixSizeX( MOCK_DATA_MATRIX_SIZE );
	mock_img.setMatrixSizeY( MOCK_DATA_MATRIX_SIZE );
	mock_img.setMatrixSizeZ( MOCK_DATA_MATRIX_SIZE );
	mock_img.setFieldOfView( MOCK_FOV, MOCK_FOV, MOCK_FOV );
	mock_img.setNumberOfChannels ( 1 );
	mock_img.setContrast( 1 );

	return mock_img;
}

ISMRMRD::Image< float > aux_test::get_mock_ismrmrd_image_with_gradients( void )
{
	size_t const Nx = MOCK_DATA_MATRIX_SIZE/2;
	size_t const Ny = MOCK_DATA_MATRIX_SIZE/2;
	size_t const Nz = MOCK_DATA_MATRIX_SIZE;
	size_t const Nc = MOCK_DATA_NUM_CHANNELS;


	ISMRMRD::Image< float > mock_img(Nx, Ny, Nz, Nc);
	for(size_t c=0; c<Nc; c++)
	for(size_t z=0; z<Nz; z++)	
	for(size_t y=0; y<Ny; y++)	
	for(size_t x=0; x<Nx; x++){

		mock_img(x,y,z,c) =  (1*x + 10*y + 100*z) + c; 

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
		*(csm.begin() +i) = std::complex<float> (0,0);
		
	
	for(size_t nc=0; nc<MOCK_DATA_NUM_CHANNELS; nc++)
	{
		for(size_t nz=0; nz<Nz; nz++)
		{
			for(size_t ny=0; ny<Ny; ny++)
			{
				for(size_t nx=0; nx<Nx; nx++)
				{
					if( nc == 0 )
						csm(nx, ny, nz, nc) = std::complex<float> (1,0);	
									
					else
						csm(nx, ny, nz, nc) = std::complex<float> (0,0);	
				}
			}
		}
	}

	return csm;
}

ISMRMRD::Image<complex_float_t> aux_test::get_mock_gaussian_csm( std::vector<size_t> vol_dims, int const num_coils )
{

	int const num_non_zero_coils = pow( (double)2,  std::floor( log2( num_coils ) ) );

	std::vector<float> sensitivity_widths {vol_dims[0] /2.f, vol_dims[1] /2.f, vol_dims[2] /2.f };


	std::vector<float> x_range, y_range, z_range;

	for( size_t i=0; i<vol_dims[0]; i++)
		x_range.push_back( i-(float)vol_dims[0]/2.f );

	for( size_t i=0; i<vol_dims[1]; i++)
		y_range.push_back( i-(float)vol_dims[1]/2.f );

	for( size_t i=0; i<vol_dims[2]; i++)
		z_range.push_back( i-(float)vol_dims[2]/2.f );

	ISMRMRD::NDArray<float> x_grid(vol_dims), y_grid(vol_dims), z_grid(vol_dims);


	for( size_t nz=0; nz<vol_dims[2]; nz++)
	for( size_t ny=0; ny<vol_dims[1]; ny++)
	for( size_t nx=0; nx<vol_dims[0]; nx++)
	{
		x_grid(nx,ny,nz) = x_range[nx];
		y_grid(nx,ny,nz) = y_range[ny];
		z_grid(nx,ny,nz) = z_range[nz];
	}

	ISMRMRD::Image<complex_float_t> csm( vol_dims[0], vol_dims[1], vol_dims[2], num_coils );

	for( size_t i=0; i<csm.getNumberOfDataElements(); i++)
		*(csm.begin() +i) = std::complex<float> (0,0);


	std::vector<size_t> center_container_size{(size_t)3, (size_t)num_coils};
	ISMRMRD::NDArray< float > coil_centers( center_container_size );
	if( num_non_zero_coils == 1)
	{
		for(size_t i=0; i<csm.getNumberOfDataElements(); i++)
			*(csm.begin() + i) = std::complex<float>(1);
	}
	else if (num_non_zero_coils == 2)
	{

		coil_centers(0, 0) = x_range[ std::floor( x_range.size()/2) ];
		coil_centers(1, 0) = y_range[0];
		coil_centers(2, 0) = 0;


		coil_centers(0, 1) = x_range[ std::floor( x_range.size()/2) ];;
		coil_centers(1, 1) = y_range.back();
		coil_centers(2, 1) = 0;

	}
	else
	{
		size_t const coils_per_side = num_non_zero_coils / 4;
		float const side_distances = 1.f/( coils_per_side + 1.f);

		for(size_t j=0; j<coils_per_side;j++)
		{
			size_t const coil_access_index = j * coils_per_side;

			coil_centers(0, coil_access_index + 0) = x_range[ std::floor( x_range.size() * float(j+1) * side_distances) ];
			coil_centers(1, coil_access_index + 0) = y_range[0];
			coil_centers(2, coil_access_index + 0) = 0;


			coil_centers(0, coil_access_index + 1) = x_range[ std::floor( x_range.size() * float(j+1) * side_distances) ];
			coil_centers(1, coil_access_index + 1) = y_range.back();
			coil_centers(2, coil_access_index + 1) = 0;


			coil_centers(0, coil_access_index + 2) = x_range[0];
			coil_centers(1, coil_access_index + 2) = y_range[ std::floor( y_range.size() * float(j+1) * side_distances) ];
			coil_centers(2, coil_access_index + 2) = 0;


			coil_centers(0, coil_access_index + 3) = x_range.back();
			coil_centers(1, coil_access_index + 3) = y_range[ std::floor( y_range.size() * float(j+1) * side_distances) ];
			coil_centers(2, coil_access_index + 3) = 0;
		}
	}

	for(size_t i_coil=0; i_coil<num_non_zero_coils; i_coil++)
	{
		for( size_t nz=0; nz<vol_dims[2]; nz++)
		for( size_t ny=0; ny<vol_dims[1]; ny++)
		for( size_t nx=0; nx<vol_dims[0]; nx++)
		{
			float const dist_square_to_coil_center = (x_grid(nx,ny,nz) - coil_centers(0, i_coil)) * (x_grid(nx,ny,nz) - coil_centers(0, i_coil))/ (sensitivity_widths[0]*sensitivity_widths[0]) +
											  (y_grid(nx,ny,nz) - coil_centers(1, i_coil)) * (y_grid(nx,ny,nz) - coil_centers(1, i_coil))/ (sensitivity_widths[1]*sensitivity_widths[1]) + 
											  (z_grid(nx,ny,nz) - coil_centers(2, i_coil)) * (z_grid(nx,ny,nz) - coil_centers(2, i_coil))/ (sensitivity_widths[2]*sensitivity_widths[2]);

			csm(nx, ny, nz, i_coil) = std::complex<float>( exp( - dist_square_to_coil_center ), 0);											  
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

ISMRMRD::AcquisitionHeader aux_test::get_mock_acquisition_header( void )
{
	ISMRMRD::AcquisitionHeader acq_hdr;
	acq_hdr.acquisition_time_stamp = 0;
	acq_hdr.number_of_samples = MOCK_DATA_RO_OVERSAMPLING * MOCK_DATA_MATRIX_SIZE;
	acq_hdr.available_channels = MOCK_DATA_NUM_CHANNELS;
	acq_hdr.center_sample = MOCK_DATA_MATRIX_SIZE/2 - 1;

	return acq_hdr;

}

AcquisitionsVector aux_test::get_mock_acquisition_vector ( ISMRMRD::IsmrmrdHeader hdr )
{

	using namespace ISMRMRD;
	typedef unsigned short unshort;


	std::ostringstream out;
	serialize(hdr, out);
	AcquisitionsVector acq_vec(out.str());

	std::vector< Encoding > encodings = hdr.encoding;

	size_t const num_scans = encodings.size();

	for( size_t iacq=0; iacq<num_scans; iacq++)
	{	

		AcquisitionSystemInformation acq_sys_info = hdr.acquisitionSystemInformation();
		unshort const NChannels = acq_sys_info.receiverChannels();

		SequenceParameters seq_par = hdr.sequenceParameters();
		std::vector<float> TE = seq_par.TE();
		std::vector<float> TR = seq_par.TR();

		Encoding enc = encodings[iacq];
		EncodingSpace enc_spac = enc.encodedSpace;

		MatrixSize size_enc_space = enc_spac.matrixSize;

		unshort const NContrasts = TE.size();
		unshort const NSamples = size_enc_space.x;
		unshort const NPhases = size_enc_space.y;
		unshort const NSlices = size_enc_space.z;
	
		for( unshort iSlice=0; iSlice<NSlices; iSlice++ )
		{
			Acquisition acq(NSamples, NChannels);
			
			if ( iSlice == 0)
				{
					acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP2);	
					
				}
			else if( iSlice == NSlices-1 )
				{
					acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2);	
				}	

			for( unshort iPhase=0; iPhase<NPhases; iPhase++ )
			{
				if ( iPhase == 0)
				{
					acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_ENCODE_STEP1);	
					acq.setFlag(ISMRMRD::ISMRMRD_ACQ_FIRST_IN_SLICE);

				}
				else if( iPhase == NPhases-1 )
				{
					acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP1);	
					acq.setFlag(ISMRMRD::ISMRMRD_ACQ_LAST_IN_SLICE);
				}	

				size_t const num_line = iSlice*NPhases + iPhase;
				AcquisitionHeader acq_hdr = acq.getHead();
				acq_hdr.idx.kspace_encode_step_1 = iPhase;
				acq_hdr.idx.kspace_encode_step_2 = iSlice;

				for( unshort iContrast=0; iContrast<NContrasts; iContrast++)
				{
					acq_hdr.scan_counter = num_line*NContrasts + iContrast;
					acq_hdr.acquisition_time_stamp = num_line*TR[0] + TE[iContrast]; 	
					
					acq_hdr.idx.contrast = 	iContrast;

					acq.setHead( acq_hdr );	
					acq_vec.append_acquisition(acq);
										
				}
			}
		}
		

		//ISMRMRD::Acquisition acq(MOCK_DATA_MATRIX_SIZE, MOCK_DATA_NUM_CHANNELS, 0);	



	}

	return acq_vec;
}

sirf::RPETrajectoryContainer aux_test::get_mock_radial_trajectory(size_t const NRad, size_t const NAng)
{
	std::vector<size_t> traj_dims;
	traj_dims.push_back(NRad);
	traj_dims.push_back(NAng);
	traj_dims.push_back((size_t)2);
	for(int i=0; i<4; i++)
		traj_dims.push_back((size_t)1);
	
	TrajVessel mock_traj( traj_dims );

	for( size_t nr=0; nr<NRad; nr++)
		for (size_t na = 0; na < NAng; na++)
		{
		    int const r_pos = nr - NRad /2;
			float const ang_pos = na*M_PI/ NAng;
			
			float const nx = r_pos * cos( ang_pos )/ (NRad);
			float const ny = r_pos * sin( ang_pos )/ (NRad);

			mock_traj(nr, na, 0) = nx;
			mock_traj(nr, na, 1) = ny;
		}

	sirf::RPETrajectoryContainer rpe_traj;
	rpe_traj.set_trajectory( mock_traj);

	return rpe_traj;
}



SignalContainer aux_test::get_mock_motion_signal()
{
	SignalContainer signal;
	float const dt = 0.1;

	std::pair<TimeAxisType, SignalAxisType> signal_point;

	for( size_t i_time=0; i_time<MOCK_NUM_SIG_PTS; i_time++)
	{
		signal_point.first = dt * i_time;
		signal_point.second = float(i_time)/float(MOCK_NUM_SIG_PTS);
		signal.push_back( signal_point);
	}

	return signal;
}


SignalContainer aux_test::get_mock_sinus_signal( AcquisitionsVector &acq_vec, TimeAxisType const period_duration_ms)
{

	unsigned const num_sampling_points = acq_vec.number();
	std::vector< TimeAxisType > all_time_points;

	for(size_t ia=0; ia<num_sampling_points; ia++)
	{
		auto sptr_acq = acq_vec.get_acquisition_sptr(ia);
		all_time_points.push_back(sptr_acq->getHead().acquisition_time_stamp); 	
	}

	auto minmax_it = std::minmax_element(std::begin(all_time_points), std::end(all_time_points));
	
	TimeAxisType t0 = *(minmax_it.first);

	SignalContainer signal;

	for( unsigned i=0; i<num_sampling_points; i++)
	{
		std::pair<TimeAxisType, SignalAxisType> signal_point;

		signal_point.first = all_time_points[i];
		signal_point.second = (1 - cos(2*M_PI / period_duration_ms * (all_time_points[i]-t0)))/2;
		signal.push_back(signal_point);
	}

	return signal;
}

SignalContainer aux_test::get_mock_sawtooth_signal( AcquisitionsVector acq_vec, TimeAxisType const period_duration_ms)
{
	
	auto sptr_acq = acq_vec.get_acquisition_sptr(0);
	TimeAxisType t_0 = sptr_acq->getHead().acquisition_time_stamp;
			

	sptr_acq = acq_vec.get_acquisition_sptr(acq_vec.items()-1);
	TimeAxisType t_fin = sptr_acq->getHead().acquisition_time_stamp;


	unsigned const num_sampling_points = acq_vec.number();

	TimeAxisType dt = float(t_fin - t_0)/ float(num_sampling_points);

	SignalContainer signal;

	for( unsigned i=0; i<num_sampling_points; i++)
	{
		std::pair<TimeAxisType, SignalAxisType> signal_point;

		signal_point.first = t_0 + i * dt;

		signal_point.second = fmod(t_0 + i * dt, period_duration_ms) / period_duration_ms;
		signal.push_back(signal_point);
	}

	return signal;	
}

SignalContainer aux_test::get_generic_contrast_inflow_signal( sirf::AcquisitionsVector &acq_vec)
{
	ISMRMRD::Acquisition acq;
	
	acq_vec.get_acquisition(0, acq);
	TimeAxisType t_0 = acq.getHead().acquisition_time_stamp;
			
	acq_vec.get_acquisition(acq_vec.items()-1, acq);
	TimeAxisType t_fin = acq.getHead().acquisition_time_stamp;

	SignalContainer signal;

	std::pair<TimeAxisType, SignalAxisType> zero_signal_point, one_signal_point;

	zero_signal_point.first =  t_0;
	zero_signal_point.second = SignalAxisType(0); 

	signal.push_back(zero_signal_point);

	one_signal_point.first = t_fin;
	one_signal_point.second = SignalAxisType(1);

	signal.push_back(one_signal_point);

	return signal;	
}

SignalContainer aux_test::get_generic_contrast_in_and_outflow_signal( sirf::AcquisitionsVector &acq_vec )
{
	ISMRMRD::Acquisition acq;
	
	acq_vec.get_acquisition(0, acq);
	TimeAxisType t_0 = acq.getHead().acquisition_time_stamp;
			
	acq_vec.get_acquisition(acq_vec.items()-1, acq);
	TimeAxisType t_fin = acq.getHead().acquisition_time_stamp;

	TimeAxisType t_half = (t_fin + t_0)/ (TimeAxisType)2;

	SignalContainer signal;

	std::pair<TimeAxisType, SignalAxisType> zero_signal_point, one_signal_point, final_signal_point;

	zero_signal_point.first =  t_0;
	zero_signal_point.second = SignalAxisType(0); 

	signal.push_back(zero_signal_point);

	one_signal_point.first = t_half;
	one_signal_point.second = SignalAxisType(1);

	signal.push_back(one_signal_point);

	final_signal_point.first = t_fin;
	final_signal_point.second = SignalAxisType(0);

	signal.push_back(final_signal_point);

	return signal;	
}


SignalContainer aux_test::get_generic_respiratory_signal( sirf::AcquisitionsVector &acq_vec)
{
	TimeAxisType const period_duartion_ms = 3300;
	return aux_test::get_mock_sinus_signal(acq_vec, period_duartion_ms );
}


SignalContainer aux_test::get_generic_cardiac_signal( sirf::AcquisitionsVector &acq_vec)
{
	TimeAxisType const period_duartion_ms = 900;
	return aux_test::get_mock_sawtooth_signal(acq_vec, period_duartion_ms );
}




