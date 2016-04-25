import math
import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *
from pGadgets import *

##def edge_weight(u):
##    shape = u.shape
##    nx = shape[0]
##    ny = shape[1]
##    mx = 2*nx - 1
##    my = 2*ny - 1
##    w = numpy.ndarray((mx, my), dtype = numpy.float32)
##    pygadgetron.find_edges(nx, ny, u.ctypes.data, w.ctypes.data)
##    return w
##
##def smoothen(u, w):
##    shape = u.shape
##    nx = shape[0]
##    ny = shape[1]
####    v = numpy.ndarray((nx, ny), dtype = numpy.float64)
##    pygadgetron.smoothen(nx, ny, u.ctypes.data, w.ctypes.data)
##    #return v

try:
    # acquisitions will be read from this HDF file
    #file = str(input('raw data file: '))
    file = 'testdata.h5'
    input_data = MR_Acquisitions(file)

    processed_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    # extract real images from complex
    images = MR_extract_real_images(complex_images)

    pp = NoiseFilter()

    # plot obtained images
    ni = 1 #images.number()
    for i in range(ni):
        data = images.image_as_array(i)
        pylab.figure(i)
        pylab.imshow(data[0,0,:,:])
        pylab.show()
        u = data[0,0,:,:]
        pp.filter(u)
##        for iter in range(10):
##            w = edge_weight(u)
##            smoothen(u, w)
        pylab.figure(1)
        pylab.imshow(u)
##        pylab.figure(2)
##        pylab.imshow(w)
##        u = data[0,0,:,:]
##        for iter in range(100):
##            u = smoothen(u, w)
##        pylab.figure(3)
##        pylab.imshow(u)
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)

