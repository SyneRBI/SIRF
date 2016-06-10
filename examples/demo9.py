'''
GRAPPA reconstruction with the steepest descent step: illustrates
the use of Acquisition Model projections
'''

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

try:
    # acquisitions will be read from this HDF file
    file = str(input('raw data file (with apostrophys in Python2.*): '))
    input_data = MR_Acquisitions(file)

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoGadget', \
         'RemoveROOversamplingGadget']
    acq_proc = AcquisitionsProcessor(prep_gadgets)
    print('---\n pre-processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)
    pp_norm = preprocessed_data.norm()

    # compute coil sensitivity maps
    csms = MR_CoilSensitivityMaps()
    print('---\n sorting acquisitions...')
    preprocessed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(preprocessed_data)

    # perform reconstruction
    recon = MR_BasicGRAPPAReconstruction()
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    output = recon.get_output()
    # for undersampled acquisition data GRAPPA computes Gfactor images
    # in addition to reconstructed ones
    complex_images = output.select(2)

    # create acquisition model based on the acquisition parameters
    # stored in preprocessed_data and image parameters stored in complex_images
    am = MR_AcquisitionModel(preprocessed_data, complex_images)
    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to simulate acquisitions
    fwd_data = am.forward(complex_images)
    fwd_norm = fwd_data.norm()

    # compute the difference between real and simulated acquisitions
    diff = fwd_data - preprocessed_data * (fwd_norm/pp_norm)
    rr = diff.norm()/fwd_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)

    # try to improve the reconstruction by the steepest descent step
    g = am.backward(diff)
    w = am.forward(g)
    alpha = (g*g)/(w*w)
    r_complex_imgs = complex_images - g*alpha

    # get real-valued reconstructed and refined images
    print('---\n processing images...')
    images = MR_extract_real_images(complex_images)
    r_imgs = MR_extract_real_images(r_complex_imgs)
    nz = images.number()
    print('%d images reconstructed.' % nz)

    # plot images
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
        rdata = r_imgs.image_as_array(i)
        pylab.figure(i + 1)
        pylab.title('GRAPPA image')
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + 1))
        pylab.figure(i + nz + 1)
        pylab.title('refined image')
        pylab.imshow(rdata[0,0,:,:])
        print('Close Figure %d window to continue...' % (i + nz + 1))
        pylab.show()

    print('done')

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
