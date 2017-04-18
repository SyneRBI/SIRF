'''
Undersampled data tests
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
##
## Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##       http://www.apache.org/licenses/LICENSE-2.0
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

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
    if na != 160:
        print('??? wrong number of acquisitions')
        return 160*2 + 4

    ny, nc, nx = acqs.dimensions()
    if (nx, ny, nc) != (512, 160, 8):
        print('??? wrong slice dimensions')
        return 160*2 + 3

    flags_failed = 0
    encode_steps_failed = 0
    for i in range(160):
        acq = acqs.acquisition(i)
        flags = 0
        if i == 0:
            flags = 64
        elif i == 159:
            flags = 128
        elif i > 47 and i < 112:
            if (i - 48)%2 == 0:
                flags = 1048576
            else:
                flags = 524288
        if i < 48:
            encode_step = 2*i
        elif i < 112:
            encode_step = i + 48
        else:
            encode_step = (i - 32)*2            
            
        if acq.flags() != flags:
            flags_failed += 1
        if acq.idx_kspace_encode_step_1() != encode_step:
            encode_steps_failed += 1

    if flags_failed > 0:
        print('??? flags failed: %d' % flags_failed)
    if encode_steps_failed > 0:
        print('??? encode_steps failed: %d' % encode_steps_failed)
    failed = flags_failed + encode_steps_failed
    if failed == 0:
        print('+++ all acquisitions tests passed')
    return failed

def main():

    failed = 0
    eps = 1e-5

    data_path = mr_data_path()
    input_data = AcquisitionData\
        (data_path + '/simulated_MR_2D_cartesian_Grappa2.h5')
    failed += acquisitions_tests_failed(input_data)

    input_data_norm = input_data.norm()
    print('---\n acquisition data norm: %e' % input_data_norm)
    failed += test_failed(1, 147.9937, input_data_norm, 0, eps)

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)

    processed_data_norm = processed_data.norm()
    print('---\n processed acquisition data norm: %e' % processed_data_norm)
    failed += test_failed(2, 142.3433, processed_data_norm, 0, eps)

    recon = CartesianGRAPPAReconstructor()
    recon.compute_gfactors(False)
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    complex_images_norm = complex_images.norm()
    print('---\n reconstructed images norm: %e' % complex_images_norm)
    failed += test_failed(3, 996.5304, complex_images_norm, 0, eps)

    csms = CoilSensitivityData()

    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n computing coil images...')
    cis = CoilImageData()
    cis.calculate(processed_data)
    print('---\n computing sensitivity maps...')
    csms.calculate(cis)

    am = AcquisitionModel(processed_data, complex_images)

    am.set_coil_sensitivity_maps(csms)

    fwd_acqs = am.forward(complex_images)

    fwd_acqs_norm = fwd_acqs.norm()
    print('---\n reconstructed images forward projection norm %e'\
          % fwd_acqs_norm)
    failed += test_failed(4, 971.9004, fwd_acqs_norm, 0, eps)

    acqs_diff = fwd_acqs - processed_data
    rr = acqs_diff.norm()/fwd_acqs_norm
    print('---\n reconstruction residual norm (rel): %e' % rr)
    failed += test_failed(5, 0.8598439, rr, eps, 0)

    bwd_images = am.backward(processed_data)
    imgs_diff = bwd_images - complex_images
    rd = imgs_diff.norm()/complex_images.norm()
    print(\
        '---\n difference between reconstructed and back-projected images: %e'\
        % rd)
    failed += test_failed(6, 0.8665925, rd, eps, 0)

    xFy = processed_data * fwd_acqs
    Bxy = bwd_images * complex_images
    print('---\n (x, F y) = (%e, %e)' % (xFy.real, xFy.imag))
    print('= (B x, y) = (%e, %e)' % (Bxy.real, Bxy.imag))
    failed += test_failed(7, xFy.real, Bxy.real, 0, eps)

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

