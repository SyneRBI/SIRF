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

#include "sirf/Reg/NiftyResample.h"
#include "sirf/Reg/NiftiImageData3DDisplacement.h"
#include "sirf/common/GeometricalInfo.h"

#include "sirf/cDynamicSimulation/auxiliary_input_output.h"

using namespace sirf;

void DynamicSimulationDeformer::deform_contrast_generator(MRContrastGenerator& mr_cont_gen, std::vector<sirf::NiftiImageData3DDeformation<float> > &vec_displacement_fields)
{

	GadgetronImagesVector& img_data = mr_cont_gen.get_contrast_filled_volumes();
	 
	NiftyResampler<float> resampler;
    resampler.set_interpolation_type_to_cubic_spline();

	std::shared_ptr< sirf::GadgetronImageData > sptr_img_to_deform = std::move(img_data.clone());

	resampler.set_reference_image(sptr_img_to_deform);
	resampler.set_floating_image (sptr_img_to_deform);

	for( size_t i_trafo=0; i_trafo<vec_displacement_fields.size(); i_trafo++)
	{
        auto disp_trafo = std::make_shared<NiftiImageData3DDeformation<float> >( vec_displacement_fields[i_trafo] );
		resampler.add_transformation(disp_trafo);
	}

	resampler.process();

	const std::shared_ptr<const sirf::ImageData>  sptr_deformed_img = resampler.get_output_sptr();
		
	sptr_deformed_img->copy(sptr_deformed_img->begin(),
							img_data.begin(), 
							img_data.end());

	std::vector< NiftiImageData3DDeformation<float> > empty_vec_to_free_memory;
	vec_displacement_fields.swap( empty_vec_to_free_memory );
}

void DynamicSimulationDeformer::deform_contrast_generator(PETContrastGenerator& pet_cont_gen, std::vector<NiftiImageData3DDeformation<float> >& vec_displacement_fields)
{
	std::vector< STIRImageData >&  vect_img_data = pet_cont_gen.get_contrast_filled_volumes();

	if( vect_img_data.size() != 2)
		throw std::runtime_error(" Please call map_tissue before the deformation of the contrast generator. You need both activity and attenaution in the correct motion state.");
	
	for(size_t i_cont=0; i_cont<vect_img_data.size(); i_cont++)
	{
		STIRImageData &curr_img = vect_img_data[i_cont];
		deform_pet_image( curr_img, vec_displacement_fields );
	}
}

void DynamicSimulationDeformer::deform_pet_image(STIRImageData& img, std::vector<NiftiImageData3DDeformation<float> >& vec_deformation_fields)
{
    std::cout << "printing vox geo inf for stir image" << std::endl;
	print_io::print_voxelized_geometrical_info( img );            

	std::shared_ptr<NiftiImageData3D<float> > sptr_img_to_deform =
            std::make_shared<NiftiImageData3D<float> >(img);
    std::cout << "#####################################################################" << std::endl;        

    std::cout << "printing vox geo inf for niftiimage from stir image" << std::endl;
	print_io::print_voxelized_geometrical_info( *sptr_img_to_deform );            

	std::cout << " \n printing vox geo inf for deformation field" << std::endl;
	print_io::print_voxelized_geometrical_info( vec_deformation_fields[0] );



	auto geom_stir_img = img.get_geom_info_sptr();

	std::vector<NiftiImageData3DDeformation<float> > deformation_fields_with_stir_geometry;
	for(size_t i=0; i<vec_deformation_fields.size(); i++)
	{
		auto sptr_nifti = vec_deformation_fields[i].get_raw_nifti_sptr();
		deformation_fields_with_stir_geometry.push_back( NiftiImageData3DDeformation<float>((float*)sptr_nifti->data, *geom_stir_img));
	}

    NiftyResampler<float> resampler;

    resampler.set_interpolation_type_to_cubic_spline();

    resampler.set_reference_image(sptr_img_to_deform);
	resampler.set_floating_image(sptr_img_to_deform);

	for( size_t i_disp=0; i_disp<deformation_fields_with_stir_geometry.size(); i_disp++)
	{
    	 std::shared_ptr<NiftiImageData3DDeformation<float> > disp_trafo =
    	 std::make_shared<NiftiImageData3DDeformation<float> >( deformation_fields_with_stir_geometry[i_disp] );
		 resampler.add_transformation(disp_trafo);
	}

	resampler.process();

}