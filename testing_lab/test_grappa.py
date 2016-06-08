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
    print('---\n processed acquisitions norm: %e' % preprocessed_data.norm())

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
    #am = MR_AcquisitionModel(preprocessed_data, complex_images)
    am = MR_AcquisitionModel(input_data, complex_images)
    am.set_coil_sensitivity_maps(csms)
    # use the acquisition model (forward projection) to produce 'acquisitions'
    fwd_data = am.forward(complex_images)
    fwd_norm = fwd_data.norm()
    print('---\n their forward projection norm %e' % fwd_norm)

    # compute the difference between real and modelled acquisitions
    #diff = fwd_data - preprocessed_data
    diff = fwd_data - input_data * (fwd_norm/input_norm)
    #c = fwd_data.norm()/preprocessed_data.norm()
    #diff = fwd_data - preprocessed_data * c
    rr = diff.norm()/fwd_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)

    # post-process reconstructed images
    print('processing images...')
    images = MR_extract_real_images(complex_images)
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
        gdata = gfactors.image_as_array(i)
        pylab.figure(i + 1)
        pylab.title('image')
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + 1))
        pylab.figure(i + nz + 1)
        pylab.title('G factor')
        pylab.imshow(gdata[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + nz + 1))
        pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
