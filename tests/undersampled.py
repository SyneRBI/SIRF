'''
Undersampled data tests
'''

import argparse
import os
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

def test_failed(ntest, expected, actual, abstol, reltol):
    if abs(expected - actual) < abstol + reltol*expected:
        print('+++ test %d passed' % ntest)
        return 0
    else:
        print('+++ test %d failed' % ntest)
        return 1

def main():

    failed = 0

    input_data = MR_Acquisitions('testdata_a2.h5')
    #failed += acquisitions_tests_failed(input_data)

    input_data_norm = input_data.norm()
    print('---\n acquisition data norm: %e' % input_data_norm)
    failed += test_failed(1, 147.9937, input_data_norm, 0, 1e-5)

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)

    processed_data_norm = processed_data.norm()
    print('---\n processed acquisition data norm: %e' % processed_data_norm)
    failed += test_failed(2, 142.3433, processed_data_norm, 0, 1e-5)

    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    complex_images_norm = complex_images.norm()
    print('---\n reconstructed images norm: %e' % complex_images_norm)
    failed += test_failed(3, 114.1564, complex_images_norm, 0, 1e-5)

    csms = MR_CoilSensitivityMaps()

    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n computing sensitivity maps...')
    csms.calculate(processed_data)

    am = MR_AcquisitionModel(processed_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    fwd_acqs = am.forward(complex_images)

    fwd_acqs_norm = fwd_acqs.norm()
    print('---\n reconstructed images forward projection norm %e'\
          % fwd_acqs_norm)
    failed += test_failed(4, 112.8888, fwd_acqs_norm, 0, 1e-5)

    acqs_diff = fwd_acqs - processed_data
    rr = acqs_diff.norm()/fwd_acqs_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)
    failed += test_failed(5, 0.7259499, rr, 1e-6, 0)

    bwd_images = am.backward(processed_data)
    imgs_diff = bwd_images - complex_images
    rd = imgs_diff.norm()/complex_images.norm()
    print(\
        '---\n difference between reconstructed and back-projected images: %e'\
        % rd)
    failed += test_failed(6, 0.6749271, rd, 1e-6, 0)

    xFy = processed_data * fwd_acqs
    Bxy = bwd_images * complex_images
    print('---\n (x, F y) = (%e, %e)' % (xFy.real, xFy.imag))
    print('= (B x, y) = (%e, %e)' % (Bxy.real, Bxy.imag))
    failed += test_failed(7, xFy.real, Bxy.real, 0, 1e-6)
##    failed += test_failed(8, 0, xFy.imag/xFy.real, 1e-6, 0)
##    failed += test_failed(9, 0, Bxy.imag/Bxy.real, 1e-6, 0)

    images = complex_images.real()
    data = images.image_as_array(0)
    #print(data[0, 0, 142, 130])
    failed += test_failed(10, 0.5952597, data[0, 0, 142, 130], 0, 1e-6)

    if failed == 0:
        print('all tests passed')
    else:
        print('%d tests failed' % failed)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
