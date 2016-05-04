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

def rel_diff(u, v):
    a = (u * v) / (v * v)
    w = u - a * v
    return w.norm()/u.norm()

def ndarray_diff(data, alt_data):
    shape = data.shape
    nx = shape[1]
    ny = shape[0]
    max_data = 0
    max_diff = 0
    xmax = -1
    ymax = -1
    for ix in range(nx):
        for iy in range(ny):
            if data[iy,ix] > max_data:
                max_data = data[iy,ix]
            diff = abs(data[iy,ix] - alt_data[iy,ix])
            if diff > max_diff:
                xmax = ix
                ymax = iy
                max_diff = diff
##    for ix in range(xmax - 1, xmax + 2):
##        for iy in range(ymax - 1, ymax + 2):
##            print(ix, iy, data[iy,ix]/max_data, alt_data[iy,ix]/max_data)
    return max_diff/max_data

try:
    # acquisitions will be read from this HDF file
    file = str(input('raw data file: '))
    input_data = MR_Acquisitions(file)

    print('---\n acquisition data norm: %e' % input_data.norm())

    processed_data = MR_remove_x_oversampling(input_data)

    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    print('---\n reconstructed images norm: %e' % complex_images.norm())

    csms = MR_CoilSensitivityMaps()

    print('---\n sorting acquisitions...')
    input_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(input_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in interim_images
    am = MR_AcquisitionModel(input_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce 'acquisitions'
    acqs = am.forward(complex_images)
    imgs = am.backward(acqs)
    im_diff = imgs - complex_images
    print('---\n 0.0 = %e' % (im_diff.norm()/imgs.norm()))

    alt_complex_images = am.backward(input_data)
##    print('---\n reconstructed images norm: %e' % alt_complex_images.norm())
##    print('---\n relative difference in complex images: %e' \
##          % rel_diff(alt_complex_images, complex_images))

    b = complex_images.norm()/alt_complex_images.norm()
    alt_complex_images = b * alt_complex_images
##    print('---\n reconstructed images norm: %e' % alt_complex_images.norm())

    # extract real images from complex
    images = MR_extract_real_images(complex_images)
    alt_images = MR_extract_real_images(alt_complex_images)

    # plot obtained images
    for i in range(images.number()):
        data = images.image_as_array(i)
        alt_data = alt_images.image_as_array(i)
        rel_diff = ndarray_diff(data[0,0,:,:], alt_data[0,0,:,:])
        print('relative maximal difference: %e' % rel_diff)
        pylab.figure(i + 1)
        pylab.imshow(alt_data[0,0,:,:])
##        pylab.figure(i + 101)
##        pylab.imshow(alt_data[0,0,:,:])
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
