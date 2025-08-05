# -*- coding: utf-8 -*-
"""sirf.STIR tests for ScatterEstimation
v{version}

Usage:
  tests_scatter [--help | options]

Options:
  -r, --record   record the measurements rather than check them [actually not used for current test]
  -v, --verbose  report each test status

{author}

{licence}
"""
import numpy
import sirf.STIR as pet
import os
from sirf.Utilities import runner, __license__
__version__ = "0.1.0"
__author__ = "Kris Thielemans"


def test_main(rec=False, verb=False, throw=True, no_ret_val=True):
    _ = pet.MessageRedirector()

    #data_path = pet.examples_data_path('PET')
    #raw_data_file = pet.existing_filepath(data_path,'Utahscat600k_ca_seg4.hs')
    data_path = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'examples', 'parameter_files')
    print(data_path)
    acq_template_filename = pet.existing_filepath(data_path, 'scatter_template.hs')
    #image_file = pet.existing_filepath(data_path, 'test_image_PM_QP_6.hv')

    #acq_data = pet.AcquisitionData(raw_data_file)
    acq_template = pet.AcquisitionData(acq_template_filename)
    act_image = pet.ImageData(acq_template)
    atten_image = act_image.get_uniform_copy(0)
    image_size = atten_image.dimensions()
    voxel_size = atten_image.voxel_sizes()

    # create a cylindrical water shape
    water_cyl = pet.EllipticCylinder()
    water_cyl.set_length(image_size[0]*voxel_size[0])
    water_cyl.set_radii((100, 100))
    water_cyl.set_origin((image_size[0]*voxel_size[0]*0.5, 0, 0))

    # add the shape to the image
    atten_image.add_shape(water_cyl, scale = 9.687E-02) # use mu of water
    act_image.add_shape(water_cyl, scale = 20)
    # Need to create ACFs
    # create acquisition model
    am = pet.AcquisitionModelUsingRayTracingMatrix()
    am.set_up(acq_template, atten_image)
    # create acquisition sensitivity model from attenuation image
    asm = pet.AcquisitionSensitivityModel(atten_image, am)
    asm.set_up(acq_template)
    am.set_acquisition_sensitivity(asm)
    # apply attenuation to the uniform acquisition data to obtain ACFs
    ACFs = acq_template.get_uniform_copy(1)
    asm.normalise(ACFs)

    sss = pet.SingleScatterSimulator()
    sss.set_attenuation_image(atten_image)
    sss.set_up(acq_template, act_image)
    scatter_data = sss.forward(act_image)

    acq_model = pet.AcquisitionModelUsingRayTracingMatrix()
    acq_model.set_acquisition_sensitivity(asm)
    acq_model.set_up(acq_template, act_image)
    unscattered_data = acq_model.forward(act_image)

    acq_data = unscattered_data + scatter_data

    # I get around 21% scatter fraction for this data
    scatter_fraction = scatter_data.norm()/acq_data.norm()
    if scatter_fraction < .18 or scatter_fraction > .25:
        scatter_data.write("out_scatter_data.hs")
        unscattered_data.write("out_unscattered.hs")
        assert False, f"Scatter fraction ({scatter_fraction}) is out of range (should be around .2 for this data)"

    scat_est = pet.ScatterEstimator()
    scat_est.set_input(acq_data)
    scat_est.set_attenuation_image(atten_image)
    scat_est.set_attenuation_correction_factors(ACFs)
    scat_est.set_asm(pet.AcquisitionSensitivityModel(acq_data.get_uniform_copy(1)))
    scat_est.set_randoms(acq_data.get_uniform_copy(0))
    scat_est.set_OSEM_num_subsets(4)
    assert scat_est.get_OSEM_num_subsets() == 4
    scat_est.set_OSEM_num_subiterations(3)
    assert scat_est.get_OSEM_num_subiterations() == 3
    scat_est.set_num_iterations(5)
    assert scat_est.get_num_iterations() == 5
    scat_est.set_up()
    scat_est.process()
    scatter_estimate = scat_est.get_output()

    rel_err = (scatter_data - scatter_estimate).norm() / scatter_estimate.norm()
    # I get around 10% error, due to randomness etc, so set our threshold a bit larger
    if  rel_err > .2:
        scatter_estimate.write("out_scatter_estimate.hs")
        scatter_data.write("out_scatter_data.hs")
        unscattered_data.write("out_unscattered.hs")
        assert False, f"Difference between simulated and estimated scatter is too large (rel err {rel_err}). Data written to file as out*.hs"

    if no_ret_val:
        return
    return 0, 1


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__, no_ret_val=False)
