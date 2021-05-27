# -*- coding: utf-8 -*-
"""sirf.STIR tests
v{version}

Usage:
  tests_NiftyPET [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import sirf.STIR as pet
from sirf.Utilities import is_operator_adjoint, runner, __license__, examples_data_path
import numpy as np
import time
import sys
import os
import unittest

has_niftypet = hasattr(pet, 'AcquisitionModelUsingNiftyPET')

__version__ = "0.2.3"
__author__ = "Richard Brown"

pet.AcquisitionData.set_storage_scheme('file')

def get_elliptical_cylinder(radius_x, radius_y, length, origin=None):
    cyl = pet.EllipticCylinder()
    cyl.set_radius_x(radius_x)
    cyl.set_radius_y(radius_y)
    cyl.set_length(length)
    if origin is not None:
        cyl.set_origin(origin)
    return cyl

def get_image():
    im_size = (127, 320, 320)
    im_spacing = (2.03125, 2.08626, 2.08626)
    image = pet.ImageData()
    image.initialise(im_size, im_spacing)
    image.fill(0)

    cyl = get_elliptical_cylinder(200,100,1000)
    image.add_shape(cyl, scale=0.75)
    cyl = get_elliptical_cylinder(100,50,300,(20,30,10))
    image.add_shape(cyl, scale=3)

    cyl = get_elliptical_cylinder(10,150,700,(-20,50,50))
    image.add_shape(cyl, scale=1.5)
    return image

def add_noise(proj_data,noise_factor = 1):
    proj_data_arr = proj_data.as_array() / noise_factor
    # Data should be >=0 anyway, but add abs just to be safe
    proj_data_arr = np.abs(proj_data_arr)
    noisy_proj_data_arr = np.random.poisson(proj_data_arr).astype('float32');
    noisy_proj_data = proj_data.clone()
    noisy_proj_data.fill(noisy_proj_data_arr);
    return noisy_proj_data
@unittest.skipUnless(has_niftypet, "NiftyPET not installed")
def test_main(rec=False, verb=False, throw=True):

    # Set STIR verbosity to off
    original_verb = pet.get_verbosity()
    pet.set_verbosity(1)

    time.sleep(0.5)
    sys.stderr.write("Testing NiftyPET projector...")
    time.sleep(0.5)

    

    # Get image
    image = get_image()

    # Get AM
    try:
        acq_model = pet.AcquisitionModelUsingNiftyPET()
    except:
        return 1, 1
    acq_model.set_cuda_verbosity(verb)

    data_path = examples_data_path('PET')
    # raw_data_file = pet.existing_filepath(data_path, 'mMR/mMR_template_span11.hs')
    raw_data_file = os.path.join(data_path, 'mMR')
    os.chdir(raw_data_file)
    template_acq_data = pet.AcquisitionData('mMR_template_span11.hs')

    acq_model.set_up(template_acq_data, image)

    # Test operator adjointness
    if verb:
        print('testing adjointness')
    if not is_operator_adjoint(acq_model, num_tests=1, verbose=True):
        raise AssertionError('NiftyPet AcquisitionModel is not adjoint')

    # Generate test data
    simulated_acq_data = acq_model.forward(image)
    simulated_acq_data_w_noise = add_noise(simulated_acq_data,10)

    obj_fun = pet.make_Poisson_loglikelihood(template_acq_data)
    obj_fun.set_acquisition_model(acq_model)

    recon = pet.OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(1)
    recon.set_num_subiterations(1)
    recon.set_input(simulated_acq_data_w_noise)
    if verb:
        print('setting up, please wait...')
    initial_estimate = image.get_uniform_copy()
    recon.set_up(initial_estimate)

    if verb:
        print('reconstructing...')
    recon.set_current_estimate(initial_estimate)
    recon.process()
    reconstructed_im = recon.get_output()
    if not reconstructed_im:
        raise AssertionError()

    # Reset original verbose-ness
    pet.set_verbosity(original_verb)

    return 0, 1

if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
