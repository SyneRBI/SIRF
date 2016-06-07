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
    ave_data = 0
    ave_diff = 0
    for ix in range(nx):
        for iy in range(ny):
            s = abs(data[iy,ix])
            if s > max_data:
                max_data = s
            diff = abs(data[iy,ix] - alt_data[iy,ix])
            ave_data = ave_data + s*s
            ave_diff = ave_diff + diff*diff
            if diff > max_diff:
                xmax = ix
                ymax = iy
                max_diff = diff
##    for ix in range(xmax - 1, xmax + 2):
##        for iy in range(ymax - 1, ymax + 2):
##            print(ix, iy, data[iy,ix]/max_data, alt_data[iy,ix]/max_data)
    ave_data = math.sqrt(ave_data)
    ave_diff = math.sqrt(ave_diff)
    #print(ave_data, ave_diff, ave_diff/ave_data)
    return max_diff/max_data, ave_diff/ave_data

try:
    # acquisitions will be read from this HDF file
    file = str(input('raw data file: '))
    input_data = MR_Acquisitions(file)

    print('---\n acquisition data norm: %e' % input_data.norm())

    processed_data = MR_remove_x_oversampling(input_data)
    #processed_data = input_data

    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    print('---\n reconstructed images norm: %e' % complex_images.norm())

    csms = MR_CoilSensitivityMaps()

    print('---\n sorting acquisitions...')
    #processed_data.sort()
    input_data.sort()
    print('---\n computing sensitivity maps...')
    #csms.calculate(processed_data)
    csms.calculate(input_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in interim_images
    am = MR_AcquisitionModel(input_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce 'acquisitions'
    acqs = am.forward(complex_images)
    diff = acqs - input_data
    rr = diff.norm()/acqs.norm()
    #rr = rel_diff(acqs, input_data)
    print('---\n reconstruction residual norm (rel): %e' % rr)

    na = input_data.number()
    nx, ny, nc = input_data.dimensions()
    nz = na//ny

    print('Enter z-coordinate of the slice to view the acquired data for it')
    print('(a value outside the range [0 : %d) will stop this loop)' % nz)
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
##        data = abs(input_data.slice_as_array(z))
##        pdata = abs(acqs.slice_as_array(z))
        data = input_data.slice_as_array(z)
        pdata = acqs.slice_as_array(z)
        diff_max, diff_ave = ndarray_diff(data[0,:,:], pdata[0,:,:])
        print('relative maximal difference: %e' % diff_max)
        print('relative average difference: %e' % diff_ave)
##        print('Enter coil number to view the acquired data for it')
##        print('(a value outside the range [0 : %d) will stop this loop)' % nc)
##        while True:
##            s = str(input('coil: '))
##            if len(s) < 1:
##                break
##            c = int(s)
##            if c < 0 or c >= nc:
##                break
##            pylab.figure(c)
##            pylab.title('input data')
##            pylab.imshow(data[c,:,:])
##            pylab.figure(c + nc)
##            pylab.title('am data')
##            pylab.imshow(pdata[c,:,:])
##            print('Close Figures %d and %d windows to continue...'% (c, c + nc))
##            pylab.show()

    print('---\n acquisition data norm: %e' % input_data.norm())

    #imgs = am.backward(acqs)
    imgs = am.backward(input_data)
    print('---\n back projected images norm: %e' % imgs.norm())

    im_diff = imgs - complex_images
    print('---\n 0.0 = %e' % (im_diff.norm()/imgs.norm()))

    acqs = am.forward(imgs)
    diff = acqs - input_data
    rr = diff.norm()/acqs.norm()
    #rr = rel_diff(acqs, input_data)
    print('---\n reconstruction residual norm (rel): %e' % rr)

##    alt_complex_images = am.backward(input_data)
####    print('---\n reconstructed images norm: %e' % alt_complex_images.norm())
####    print('---\n relative difference in complex images: %e' \
####          % rel_diff(alt_complex_images, complex_images))
##
##    b = complex_images.norm()/alt_complex_images.norm()
##    alt_complex_images = b * alt_complex_images
####    print('---\n reconstructed images norm: %e' % alt_complex_images.norm())
##
##    # extract real images from complex
##    images = MR_extract_real_images(complex_images)
##    alt_images = MR_extract_real_images(alt_complex_images)
##
##    # plot obtained images
##    for i in range(images.number()):
##        data = images.image_as_array(i)
##        alt_data = alt_images.image_as_array(i)
##        rel_diff = ndarray_diff(data[0,0,:,:], alt_data[0,0,:,:])
##        print('relative maximal difference: %e' % rel_diff)
##        pylab.figure(i + 1)
##        pylab.imshow(alt_data[0,0,:,:])
####        pylab.figure(i + 101)
####        pylab.imshow(alt_data[0,0,:,:])
##        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
