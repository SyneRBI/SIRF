import argparse
import numpy
import os
import pylab
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
OSMAPOSL reconstruction demo with all parameters defined in the script
and user-controlled iterations
''')
args = parser.parse_args()

def main():

    # direct all information printing to a file
    info_printer = printerTo('test2i.txt', INFO_CHANNEL)
    # write all warnings in a file
    warn_printer = printerTo('test2w.txt', WARNING_CHANNEL)

    # create acquisition model
    am = PETAcquisitionModelUsingMatrix()
    am.set_matrix(RayTracingMatrix())

    # create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqModData()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(PETAcquisitionData('my_raw_data.hs'))

    num_subiterations = 2

    # create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subiterations(num_subiterations)
    recon.set_num_subsets(num_subiterations)
    recon.set_output_filename_prefix('reconstructedImage')

    # create initial image estimate
    image = PETImage()
    image.initialise((100, 100, 30), (3, 3, 3.375))
    image.fill(1.0)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # run monitored reconstruction
    slice = 20
    for iter in range(1, num_subiterations + 1):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        # plot the current image
        data = image.as_array() # 3D NumPy array
        pylab.figure(iter)
        pylab.imshow(data[slice,:,:])
        print('close Figure %d window to continue' % iter)
        pylab.show()

    # compare the reconstructed image to the exact image
    exactImage = PETImage('my_image.hv')
    x_data = exactImage.as_array()
    data = image.as_array()
    pylab.figure(1000)
    pylab.imshow(data[slice,:,:])
    pylab.figure(1001)
    pylab.imshow(x_data[slice,:,:])
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
