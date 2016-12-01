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

def my_image_data_processor(image_array):
    return image_array

def main():

    # direct all information printing and warnings to files
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)

    # define acquisition model as one whereby the geometric
    # forward projection is represented by a ray tracing matrix
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())

    # define acquisition data source
    ad = AcquisitionData('my_forward_projection.hs')

    # define initial image estimate compatible with raw data
    # and initialize it to value 1.0 at each voxel
    image = ad.create_empty_image(1.0)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood with linear model for mean
    # and using acquisition model
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    num_subiterations = 2

    # define reconstruction algorithm to be used as
    # Ordered Subsets Maximum A-Priori One Step Late
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(12)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    recon.set_current_estimate(image)
    for iter in range(num_subiterations):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update_current_estimate()
        # copy current image estimate into python array to inspect/process
        image_array = recon.get_current_estimate().as_array()
        # plot the current image
        pylab.figure(iter + 1)
        pylab.imshow(image_array[20,:,:])
        print('close Figure %d window to continue' % (iter + 1))
        pylab.show()
        # apply your image data processor
        my_image_array = my_image_data_processor(image_array)
        # fill the current image estimate with new data
        recon.get_current_estimate().fill(my_image_array)

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
