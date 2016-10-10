import argparse
import numpy
import pylab
import os
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')
import time

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
OSMAPOSL reconstruction demo with all parameters defined in the script
and user-controlled iterations
''')
args = parser.parse_args()

def main():

    # direct all information printing to a file
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    # direct all warning printing to a file
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = printerTo('stdout', ERROR_CHANNEL)

    # create matrix to be used by the acquisition model
    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # define acquisition data
    ad = AcquisitionData('my_forward_projection.hs')

    # create prior
    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.001)

    # create filter
    filter = CylindricFilter()

    # create initial image estimate
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375)
    image = Image()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
##    filter.set_strictly_less_than_radius(False)
##    filter.apply(image)
##    filter.set_strictly_less_than_radius(True)

    # create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

    num_subiterations = 6

    # create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_MAP_model('multiplicative')
    recon.set_num_subsets(12)
    recon.set_num_subiterations(num_subiterations)
    recon.set_save_interval(num_subiterations)
    recon.set_inter_iteration_filter_interval(1)
    recon.set_inter_iteration_filter(filter)
    recon.set_output_filename_prefix('reconstructedImage')

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # plot the initial image
    data = image.as_array()
    pylab.figure(1)
    pylab.imshow(data[20,:,:])
    print('Figure 1: initial image - close window to continue')
    pylab.show()

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iter in range(1, num_subiterations + 1):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        # plot the current image
        data = image.as_array()
        pylab.figure(iter + 1)
        pylab.imshow(data[20,:,:])
        print('close Figure %d window to continue' % (iter + 1))
        pylab.show()
        # image can be post-processed
        #filter.apply(image)

    # compare the reconstructed image to the expected image
    expectedImage = Image('expected_image.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image: %e' % diff)

    # compare the reconstructed image to the exact image
    exactImage = Image('my_image.hv')
    x_data = exactImage.as_array()
    data = image.as_array()

    pylab.figure(1000)
    pylab.imshow(data[20,:,:])
    pylab.figure(1001)
    pylab.imshow(x_data[20,:,:])
    print('close Figure 1001 window, then')
    print('close Figure 1000 window to continue')
    pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
