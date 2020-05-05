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
from sirf.STIR import *
from sirf.Utilities import is_operator_adjoint, runner, __license__
import numpy as np

__version__ = "0.2.3"
__author__ = "Richard Brown"


def get_elliptical_cylinder(radius_x, radius_y, length, origin=None):
    cyl = EllipticCylinder()
    cyl.set_radius_x(radius_x)
    cyl.set_radius_y(radius_y)
    cyl.set_length(length)
    if origin is not None:
        cyl.set_origin(origin)
    return cyl

def get_image():
    im_size = (127, 320, 320)
    im_spacing = (2.03125, 2.08626, 2.08626)
    image = ImageData()
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

def test_main(rec=False, verb=False, throw=True):

    # Set STIR verbosity to off
    original_verb = get_verbosity()
    set_verbosity(0)

    time.sleep(0.5)
    sys.stderr.write("Testing NiftyPET projector...")
    time.sleep(0.5)

    data_path = examples_data_path('PET')
    raw_data_file = existing_filepath(data_path, 'mMR/mMR_template_span11.hs')
    template_acq_data = AcquisitionData(raw_data_file)

    # Get image
    image = get_image()

    # Get AM
    try:
        acq_model = AcquisitionModelUsingNiftyPET()
    except:
        return 1, 1
    acq_model.set_cuda_verbosity(verb)
    acq_model.set_up(template_acq_data, image)

    # Test operator adjointness
    if verb:
        print('testing adjointness')
    if not is_operator_adjoint(acq_model, verbose = verb):
        raise AssertionError('NiftyPet AcquisitionModel is not adjoint')

    # Generate test data
    simulated_acq_data = acq_model.forward(image)
    simulated_acq_data_w_noise = add_noise(simulated_acq_data,10)

    obj_fun = make_Poisson_loglikelihood(template_acq_data)
    obj_fun.set_acquisition_model(acq_model)

    recon = OSMAPOSLReconstructor()
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
    set_verbosity(original_verb)

    return 0, 1

if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
