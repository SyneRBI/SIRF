'''
Fully sampled data tests
'''

import argparse
import os
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'
DATA_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/examples/'

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

def acquisitions_tests_failed(acqs):

    na = acqs.number()
    if na != 512:
        print('??? wrong number of acquisitions')
        return 512*3 + 4

    nx, ny, nc = acqs.slice_dimensions()
    if (nx, ny, nc) != (512, 256, 8):
        print('??? wrong slice dimensions')
        return 512*3 + 3

    flags_failed = 0
    encode_steps_failed = 0
    repetitions_failed = 0
    for i in range(512):
        acq = acqs.acquisition(i)
        flags = 0
        if i == 0 or i == 256:
            flags = 64
        elif i == 255 or i == 511:
            flags = 128
        if acq.flags() != flags:
            flags_failed += 1
        if acq.idx_kspace_encode_step_1() != i%256:
            encode_steps_failed += 1
        if acq.idx_repetition() != i//256:
            repetitions_failed += 1

    if flags_failed > 0:
        print('??? flags failed: %d' % flags_failed)
    if encode_steps_failed > 0:
        print('??? encode_steps failed: %d' % encode_steps_failed)
    if repetitions_failed > 0:
        print('??? repetitions failed: %d' % repetitions_failed)
    failed = flags_failed + encode_steps_failed + repetitions_failed
    if failed == 0:
        print('+++ all acquisitions tests passed')
    return failed

def main():

    failed = 0
    eps = 1e-4

    input_data = AcquisitionData(DATA_PATH + 'simulated_MR_2D_cartesian.h5')
    failed += acquisitions_tests_failed(input_data)

    input_data_norm = input_data.norm()
    print('---\n acquisition data norm: %e' % input_data_norm)
    failed += test_failed(1, 221.2011, input_data_norm, 0, eps)

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)

    processed_data_norm = processed_data.norm()
    print('---\n processed acquisition data norm: %e' % processed_data_norm)
    failed += test_failed(2, 209.021, processed_data_norm, 0, eps)

    recon = SimpleReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    complex_images_norm = complex_images.norm()
    print('---\n reconstructed images norm: %e' % complex_images_norm)
    failed += test_failed(3, 209.021, complex_images_norm, 0, eps)

    cis = CoilImages()

    csms = CoilSensitivityMaps()

    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n computing coil images...')
    cis.calculate(processed_data)
    print('---\n computing sensitivity maps...')
    csms.calculate(cis) #, method = 'Inati(niter = 5)')

    am = AcquisitionModel(processed_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    fwd_acqs = am.forward(complex_images)

    fwd_acqs_norm = fwd_acqs.norm()
    print('---\n reconstructed images forward projection norm %e'\
          % fwd_acqs_norm)
    failed += test_failed(4, 209.021, fwd_acqs_norm, 0, eps)

    acqs_diff = fwd_acqs - processed_data
    rr = acqs_diff.norm()/fwd_acqs_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)
    failed += test_failed(5, 0, rr, 1e-6, 0)

    bwd_images = am.backward(processed_data)
    imgs_diff = bwd_images - complex_images
    rd = imgs_diff.norm()/complex_images.norm()
    print(\
        '---\n difference between reconstructed and back-projected images: %e'\
        % rd)
    failed += test_failed(6, 0, rd, 1e-6, 0)

    xFy = processed_data * fwd_acqs
    Bxy = bwd_images * complex_images
    print('---\n (x, F y) = (%e, %e)' % (xFy.real, xFy.imag))
    print('= (B x, y) = (%e, %e)' % (Bxy.real, Bxy.imag))
    failed += test_failed(7, xFy.real, Bxy.real, 0, 1e-6)
    failed += test_failed(8, 0, xFy.imag/xFy.real, 1e-6, 0)
    failed += test_failed(9, 0, Bxy.imag/Bxy.real, 1e-6, 0)

    if failed == 0:
        print('all tests passed')
    else:
        print('%d tests failed' % failed)
    return failed

try:
    failed = main()
    print('done')
    if failed != 0:
        sys.exit(failed)

except error as err:
    # display error information
    print('??? %s' % err.value)
    sys.exit(-1)

