# -*- coding: utf-8 -*-
'''Test set 2.

Undersampled data tests

Usage:
  test1 [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
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

"""
Created on Tue Nov 21 11:23:39 2017

@author: Evgueni Ovtchinnikov
"""

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from pGadgetron import *

record = args['--record']
verbose = args['--verbose']

def main(rec = record, verb = verbose):

    test = pTest('test2.txt', rec)
    test.verbose = verb

    data_path = mr_data_path()
    input_data = AcquisitionData\
        (data_path + '/simulated_MR_2D_cartesian_Grappa2.h5')
    test.check(input_data.norm())

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    test.check(processed_data.norm())

    recon = CartesianGRAPPAReconstructor()
    recon.compute_gfactors(False)
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()
    test.check(complex_images.norm())

    csms = CoilSensitivityData()

    processed_data.sort()
    cis = CoilImageData()
    cis.calculate(processed_data)
    csms.calculate(cis)

    am = AcquisitionModel(processed_data, complex_images)
    am.set_coil_sensitivity_maps(csms)
    fwd_acqs = am.forward(complex_images)
    fwd_acqs_norm = fwd_acqs.norm()
    test.check(fwd_acqs_norm)

    acqs_diff = fwd_acqs - processed_data
    rr = acqs_diff.norm()/fwd_acqs_norm
    test.check(rr)

    bwd_images = am.backward(processed_data)
    imgs_diff = bwd_images - complex_images
    rd = imgs_diff.norm()/complex_images.norm()
    test.check(rd)
    xFy = processed_data * fwd_acqs
    Bxy = bwd_images * complex_images
    test.check(abs(xFy.real/Bxy.real - 1), abs_tol = 1e-4)

    return test.failed, test.ntest

if __name__ == '__main__':

    try:
        failed, ntest = main()
        if failed == 0:
            if not record:
                print('all tests passed')
            else:
                print('%d measurements recorded' % ntest)
            sys.exit(0)
        else:
            print('%d of the tests failed' % failed)
            sys.exit(failed)

    except error as err:
        # display error information
        print('??? %s' % err.value)
        sys.exit(-1)
