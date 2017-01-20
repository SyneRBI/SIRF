import argparse
import numpy
import os
import pylab
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
Shortest demo not using .par file
''')
args = parser.parse_args()

def main():

    ad = AcquisitionData('my_raw_data.hs')
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    image = Image(ad).fill(1.0)
    recon.set_up(image)
    recon.reconstruct(image)

    exactImage = Image('my_image.hv')
    x_data = exactImage.as_array()
    data = image.as_array()
    pylab.figure(1000)
    pylab.imshow(data[20,:,:])
    pylab.figure(1001)
    pylab.imshow(x_data[20,:,:])
    pylab.show()

try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
