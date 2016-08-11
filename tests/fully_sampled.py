'''
Fully sampled data test
'''

import argparse
import os
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions('testdata.h5')

    print('---\n acquisition data norm: %e' % input_data.norm())

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)

    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    print('---\n reconstructed images norm: %e' % complex_images.norm())

    csms = MR_CoilSensitivityMaps()

    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(processed_data)

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in complex_images
    am = MR_AcquisitionModel(processed_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce 'acquisitions'
    acqs = am.forward(complex_images)

    print('---\n reconstructed images forward projection norm %e' % acqs.norm())

    # compute the difference between real and modelled acquisitions
    diff = acqs - processed_data
    rr = diff.norm()/acqs.norm()
    print('---\n reconstruction residual norm (rel): %e' % rr)

    # compare reconstructed images with backward-projected acquisitions
    bwd_images = am.backward(processed_data)
    im_diff = bwd_images - complex_images
    print(\
        '---\n difference between reconstructed and back-projected images: %e'\
        % (im_diff.norm()/complex_images.norm()))

    # test that the backward projection is the adjoint of forward
    # on x = processed_data and y = complex_images
    xFy = processed_data * acqs
    Bxy = bwd_images * complex_images
    print('---\n (x, F y) = (%e, %e)' % (xFy.real, xFy.imag))
    print('= (B x, y) = (%e, %e)' % (Bxy.real, Bxy.imag))

    # extract real images from complex
    images = complex_images.real()
    data = images.image_as_array(0)
    print(data[0, 0, 142, 130])

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
