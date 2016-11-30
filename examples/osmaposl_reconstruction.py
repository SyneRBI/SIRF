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

    # direct all information printing and warnings to files
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)

    # create acquisition model that uses ray tracing matrix
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())

    # define acquisition data
    ad = AcquisitionData('my_forward_projection.hs')

    # create initial image estimate compatible with raw data
    image = Image(ad).fill(1.0)

    # create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    num_subiterations = 6

    # create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(12)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iter in range(1, num_subiterations + 1):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        # copy current image estimate into python array to inspect/process;
        # image.fill(image_array) fills current estimate with processed data
        image_array = image.as_array()
        # plot the current image
        pylab.figure(iter + 1)
        pylab.imshow(image_array[20,:,:])
        print('close Figure %d window to continue' % iter)
        pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
