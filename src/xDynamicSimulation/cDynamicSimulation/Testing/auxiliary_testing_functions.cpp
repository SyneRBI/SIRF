/* ================================================

Author: Johannes Mayer
Date: 2018.03.28
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */


#include "auxiliary_testing_functions.h"


#include <omp.h>
#include <sstream>

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
	pet_tissue_pars.attenuation_1_by_mm_ = 0.01;
	pet_tissue_pars.suv_ = 15;


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
	
	tiss_at_0.pet_tissue_.attenuation_1_by_mm_= 0;
	tiss_at_0.pet_tissue_.suv_= 0;

	output.first = tiss_at_0;

	tiss_at_1.name_ = "dynamic_par";
	tiss_at_1.label_ = 10000;

	tiss_at_1.mr_tissue_.spin_density_percentH2O_ = 80;
	tiss_at_1.mr_tissue_.t1_miliseconds_ = 600;
	tiss_at_1.mr_tissue_.t2_miliseconds_ = 44;
	tiss_at_1.mr_tissue_.cs_ppm_ = 0;
	
	tiss_at_1.pet_tissue_.attenuation_1_by_mm_= 0;
	tiss_at_1.pet_tissue_.suv_= 0;

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

	pet_cont_gen.set_template_image_from_file( PET_TEMPLATE_IMAGE_DATA_PATH );

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


SignalContainer aux_test::get_mock_motion_signal( AcquisitionsVector acq_vec)
{
	SignalContainer signal;
	
	std::pair<TimeAxisType, SignalAxisType> signal_point;
	
	ISMRMRD::Acquisition acq;
	
	acq_vec.get_acquisition(0, acq);
	signal_point.first = acq.getHead().acquisition_time_stamp;
	signal_point.second = 0;

	std::cout << "(" << signal_point.first << "/" << signal_point.second << ")" <<std::endl;

	signal.push_back(signal_point);
	
	acq_vec.get_acquisition(acq_vec.items()-1, acq);
	
	signal_point.first = acq.getHead().acquisition_time_stamp;
	signal_point.second = 1;

	std::cout << "(" << signal_point.first<< "/"  << signal_point.second << ")" <<std::endl;
	signal.push_back(signal_point);

	return signal;
}