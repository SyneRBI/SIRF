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
Mutator/accessor methods demo
''')
args = parser.parse_args()

def main():

    # direct all printing to a file
    printer = printerTo('stir_demo1.txt')

    # create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction('OSMAPOSL_test_PM_QP2.par')

    # check/redefine some parameters
    f = CylindricFilter(recon.get_inter_iteration_filter())
    f.set_strictly_less_than_radius(True)
    obj = recon.get_objective_function()
    prior = obj.get_prior()
    print('prior penalisation factor: %f' % prior.get_penalisation_factor())
    prior.set_penalisation_factor(0.001)

    # read an initial estimate for the reconstructed image from a file
    image = PETImage('my_image0.hv')

    # plot the initial image
    data = image.as_array()
    pylab.figure(1)
    pylab.imshow(data[20,:,:])
    print('Figure 1: initial image - close window to continue')
    pylab.show()

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # run reconstruction
    print('reconstructing, please wait...')
    start_time = time.time()
    recon.reconstruct(image)
    elapsed_time = time.time() - start_time
    print('elapsed time: %f' % elapsed_time)

    # compare the reconstructed image to the expected image
    expectedImage = PETImage('expected_image.hv')
    diff = expectedImage.diff_from(image)
    print('difference from expected image: %e' % diff)

    # compare the reconstructed image to the exact image
    exactImage = PETImage('my_image.hv')

    # let the user inspect any z-crossections of the image they want to
    data = image.as_array()
    x_data = exactImage.as_array()
    nz = data.shape[0]
    print('Enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop the loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        pylab.figure(z)
        pylab.imshow(data[z,:,:])
        pylab.figure(1000 + z)
        pylab.imshow(x_data[z,:,:])
        print('close Figure %d window, then' % (1000 + z))
        print('close Figure %d window to continue' % z)
        pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
