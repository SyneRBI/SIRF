/* ================================================

Author: Johannes Mayer
Date: 2018.08.27
Institution: Physikalisch-Technische Bundesanstalt Berlin

================================================ */

#include "sirf/cDynamicSimulation/dynsim_deformer.h"


#include <stdexcept>
#include <sstream>
#include <complex>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>

#include "sirf/cReg/NiftyResample.h"
#include "sirf/cReg/NiftiImageData3DDisplacement.h"
#include "sirf/common/GeometricalInfo.h"

using namespace sirf;

ISMRMRD::Image< float > DynamicSimulationDeformer::extract_real_part(  ISMRMRD::Image< complex_float_t >& complex_img )
{
	return DynamicSimulationDeformer::extract_complex_subpart( complex_img, true);
}


ISMRMRD::Image< float > DynamicSimulationDeformer::extract_imaginary_part(  ISMRMRD::Image< complex_float_t >& complex_img )
{
	return DynamicSimulationDeformer::extract_complex_subpart( complex_img, false);
}

ISMRMRD::Image< float > DynamicSimulationDeformer::extract_complex_subpart( ISMRMRD::Image< complex_float_t >& complex_img, bool const extract_real_part )
{

	auto cplx_img_header = complex_img.getHead();
	cplx_img_header.data_type = 5; 

	ISMRMRD::Image< float > sub_img;
	sub_img.setHead(cplx_img_header);
	sub_img.resize( sub_img.getMatrixSizeX(), sub_img.getMatrixSizeY(), sub_img.getMatrixSizeZ(), sub_img.getNumberOfChannels() );


	for( size_t i=0; i<complex_img.getNumberOfDataElements(); i++)
	{
		if (extract_real_part)
		{
			*(sub_img.begin() + i) = std::real(  *(complex_img.begin() + i ) );
		}
		else
		{
			*(sub_img.begin() + i) = std::imag(  *(complex_img.begin() + i ) );
		}
	}

	return sub_img;
}



void DynamicSimulationDeformer::deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields)
{

	std::vector< ISMRMRD::Image< complex_float_t> >&  vect_img_data = mr_cont_gen.get_contrast_filled_volumes();

	for(size_t i_cont=0; i_cont<vect_img_data.size(); i_cont++)
	{
		ISMRMRD::Image<complex_float_t> &curr_img = vect_img_data[i_cont];
		auto real_image_part =  DynamicSimulationDeformer::extract_real_part(curr_img);
		DynamicSimulationDeformer::deform_ismrmrd_image( real_image_part, vec_displacement_fields);
		
		auto imaginary_image_part =  DynamicSimulationDeformer::extract_imaginary_part(curr_img);
		DynamicSimulationDeformer::deform_ismrmrd_image( imaginary_image_part, vec_displacement_fields);			
		
		for( size_t i_vox=0; i_vox<curr_img.getNumberOfDataElements(); i_vox++)
		{
			float const voxel_real_part = *(real_image_part.begin() + i_vox);
			float const voxel_imag_part = *(imaginary_image_part.begin() + i_vox);
			
			*(curr_img.begin() + i_vox) =  std::complex<float> ( voxel_real_part, voxel_imag_part );
		}
	}
	
	std::vector< NiftiImageData3DDeformation<float> > empty_vec_to_free_memory;
	vec_displacement_fields.swap( empty_vec_to_free_memory );
}
	


void DynamicSimulationDeformer::deform_ismrmrd_image(ISMRMRD::Image< float >& img, std::vector<NiftiImageData3DDeformation<float> > &vec_deformation_fields)
{
	
    std::shared_ptr<const VoxelisedGeometricalInfo3D > sptr_geometrical_info = vec_deformation_fields[0].get_geom_info_sptr();

    std::cout << "WARNING: TAKING THE FIRST VECTORFIELDs GEOMETRY " << std::endl;

    auto sptr_img_to_deform = std::make_shared< NiftiImageData3D<float> >(img.getDataPtr(), *sptr_geometrical_info);

	NiftyResample<float> resampler;

    resampler.set_interpolation_type_to_cubic_spline();
	resampler.set_reference_image(sptr_img_to_deform);
	resampler.set_floating_image (sptr_img_to_deform);

	for( size_t i_trafo=0; i_trafo<vec_deformation_fields.size(); i_trafo++)
	{
        auto disp_trafo = std::make_shared<NiftiImageData3DDeformation<float> >( vec_deformation_fields[i_trafo] );
		resampler.add_transformation(disp_trafo);
	}

	resampler.process();

	auto deformed_img = resampler.get_output_sptr();
    auto deformed_img_as_nifti = deformed_img->get_raw_nifti_sptr();
	
	if( deformed_img_as_nifti->nvox != img.getNumberOfDataElements() )
		throw std::runtime_error("Something went wrong during the resampling. The output image and input image have different number of voxels.");

	memcpy(img.begin(), deformed_img_as_nifti->data, sizeof(float) * deformed_img_as_nifti->nvox);
}


void DynamicSimulationDeformer::deform_contrast_generator(PETContrastGenerator& mr_cont_gen, std::vector<NiftiImageData3DDeformation<float> >& vec_displacement_fields)
{
	std::vector< STIRImageData >&  vect_img_data = mr_cont_gen.get_contrast_filled_volumes();

	if( vect_img_data.size() != 2)
		throw std::runtime_error(" Please call map_tissue before the deformation of the contrast generator. You need both activity and attenaution in the correct motion state.");
	
	for(size_t i_cont=0; i_cont<vect_img_data.size(); i_cont++)
	{
		STIRImageData &curr_img = vect_img_data[i_cont];
		deform_pet_image( curr_img, vec_displacement_fields );
	}
}

void DynamicSimulationDeformer::deform_pet_image(STIRImageData& img, std::vector<NiftiImageData3DDeformation<float> >& vec_displacement_fields)
{
	std::shared_ptr<NiftiImageData3D<float> > img_to_deform =
            std::make_shared<NiftiImageData3D<float> >(img);

    NiftyResample<float> resampler;

    resampler.set_interpolation_type_to_cubic_spline();
	resampler.set_reference_image(img_to_deform);
	resampler.set_floating_image(img_to_deform);

	for( size_t i_disp=0; i_disp<vec_displacement_fields.size(); i_disp++)
	{
        std::shared_ptr<NiftiImageData3DDisplacement<float> > disp_trafo =
                std::make_shared<NiftiImageData3DDisplacement<float> >( vec_displacement_fields[i_disp] );
		resampler.add_transformation(disp_trafo);
	}

	resampler.process();

	auto deformed_img = resampler.get_output_sptr();
    img.sirf::ImageData::fill(*deformed_img);
}