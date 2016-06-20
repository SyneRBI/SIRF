'''
Upper level 3-steps GRAPPA reconstruction demo.
'''

import math
import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'
DATA_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/examples/'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    file = str(input('raw data file (with apostrophys in Python2.*): '))
    input_data = MR_Acquisitions(DATA_PATH + file)
    input_norm = input_data.norm()
    print('---\n acquisitions norm: %e' % input_norm)

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoGadget', \
         'RemoveROOversamplingGadget']
    acq_proc = AcquisitionsProcessor(prep_gadgets)
    print('pre-processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)
    pp_norm = preprocessed_data.norm()
    print('---\n processed acquisitions norm: %e' % pp_norm)

    csms = MR_CoilSensitivityMaps()
    print('---\n sorting acquisitions...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(preprocessed_data)

    recon = MR_BasicGRAPPAReconstruction()
    # connect to input data
    recon.set_input(preprocessed_data)
    # perform reconstruction
    print('---\n reconstructing...')
    recon.process()
    # get reconstructed images
    output = recon.get_output()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    complex_images = output.select(2)
    complex_gfactors = output.select(2, 1)
    print('---\n reconstructed images norm: %e' % complex_images.norm())

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in complex_images
    am = MR_AcquisitionModel(preprocessed_data, complex_images)
    #am = MR_AcquisitionModel(input_data, complex_images)
    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce 'acquisitions'
    fwd_data = am.forward(complex_images)
    fwd_norm = fwd_data.norm()
    print('---\n their forward projection norm %e' % fwd_norm)

##    #bwd_images = am.backward(preprocessed_data)
##    bwd_images = am.backward(fwd_data)
##    c = (bwd_images.norm()/complex_images.norm())
##    im_diff = bwd_images - complex_images*c
##    print('---\n 0.0 = %e' % (im_diff.norm()/bwd_images.norm()))

    # compute the difference between real and modelled acquisitions
    #diff = fwd_data - preprocessed_data
    #diff = fwd_data - input_data * (fwd_norm/input_norm)
    diff = fwd_data - preprocessed_data * (fwd_norm/pp_norm)
    rr = diff.norm()/fwd_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)

    nx, ny, nc = input_data.slice_dimensions()
    nz = complex_images.number()
    print('Enter z-coordinate of the slice to view the acquired data for it')
    print('(a value outside the range [0 : %d) will stop this loop)' % nz)
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
##        data = abs(diff.slice_as_array(z))
        data = abs(preprocessed_data.slice_as_array(z)) * (fwd_norm/pp_norm)
        pdata = abs(fwd_data.slice_as_array(z))
        pdata = data - pdata
####        data = input_data.slice_as_array(z)
####        pdata = acqs.slice_as_array(z)
##        diff_max, diff_ave = ndarray_diff(data[0,:,:], pdata[0,:,:])
##        print('relative maximal difference: %e' % diff_max)
##        print('relative average difference: %e' % diff_ave)
        print('Enter coil number to view the acquired data for it')
        print('(a value outside the range [0 : %d) will stop this loop)' % nc)
        while True:
            s = str(input('coil: '))
            if len(s) < 1:
                break
            c = int(s)
            if c < 0 or c >= nc:
                break
            pylab.figure(c)
            pylab.title('input data')
            pylab.imshow(data[c,:,:])
            pylab.colorbar();
            pylab.figure(c + nc)
            pylab.title('diff')
            #pylab.title('am data')
            pylab.imshow(pdata[c,:,:])
            print('Close Figures %d and %d windows to continue...'% (c, c + nc))
            pylab.colorbar();
            pylab.show()

    # try to improve the reconstruction by steepest descent iterations
    x = complex_images*(pp_norm/fwd_norm)
    f = am.forward(x)
    r = f - preprocessed_data
    g = am.backward(r)
    w = am.forward(g)
    alpha = (g*g)/(w*w)
    x = x - g*alpha
    imgs = MR_extract_real_images(x)
##    gamma = g*g
##    p = g
##    for iter in range(3):
##        w = am.forward(p)
##        alpha = (g*p)/(w*w)
##        x = x - p*alpha
##        imgs = MR_extract_real_images(x)
##        idata = imgs.image_as_array(0)
##        pylab.figure(iter + 1)
##        pylab.title('iterated image')
##        pylab.imshow(idata[0,0,:,:])
##        f = f - w*alpha
##        #f = am.forward(x)
##        r = f - preprocessed_data
##        g = am.backward(r)
##        print('---\n residual norm %e, gradient norm %e' % (r.norm()/pp_norm, g.norm()))
##        beta = g*g
##        p = g + p*(beta/gamma)
##        gamma = beta
    pylab.show()

    print('Enter z-coordinate of the slice to view the acquired data for it')
    print('(a value outside the range [0 : %d) will stop this loop)' % nz)
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = abs(preprocessed_data.slice_as_array(z))
        f = am.forward(x)
        fdata = abs(f.slice_as_array(z))
        pdata = data - fdata
        print('Enter coil number to view the acquired data for it')
        print('(a value outside the range [0 : %d) will stop this loop)' % nc)
        while True:
            s = str(input('coil: '))
            if len(s) < 1:
                break
            c = int(s)
            if c < 0 or c >= nc:
                break
            pylab.figure(c)
            pylab.title('input data')
            pylab.imshow(data[c,:,:])
            pylab.colorbar();
            pylab.figure(c + nc)
            pylab.title('diff')
            pylab.imshow(pdata[c,:,:])
            print('Close Figures %d and %d windows to continue...'% (c, c + nc))
            pylab.colorbar();
            pylab.show()

    # post-process reconstructed images
    print('processing images...')
    #images = MR_extract_real_images(bwd_images)
    images = MR_extract_real_images(complex_images)
##    imgs = MR_extract_real_images(x)
    gfactors = MR_extract_real_images(complex_gfactors)

    nz = images.number()
    print('%d images reconstructed.' % nz)

    # plot obtained images
    print('Enter z-coordinate of the slice to view it')
    print('(a value outside the range [0 : %d] will stop this loop)'%(nz - 1))
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        i = z
        data = images.image_as_array(i)
        idata = imgs.image_as_array(i)
        gdata = gfactors.image_as_array(i)
        pylab.figure(i + 1)
        pylab.title('GRAPPA image')
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + 1))
        pylab.figure(i + nz + 1)
##        pylab.title('G factor')
##        pylab.imshow(gdata[0,0,:,:])
        pylab.title('iterated image')
        pylab.imshow(idata[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + nz + 1))
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
