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
    info_printer = printerTo('test1i.txt', INFO_CHANNEL)
    # write all warnings in a file
    warn_printer = printerTo('test1w.txt', WARNING_CHANNEL)

    # define acquisition data
    ad = AcquisitionData('my_raw_data.hs')

    # define acquisition model
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())

    # define objective function
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    # define reconstructor
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subiterations(2)
    recon.set_num_subsets(2)

    # create initial image estimate
    image = Image(ad).fill(1.0)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # run reconstruction
    print('reconstructing, please wait...')
    recon.reconstruct(image)

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
