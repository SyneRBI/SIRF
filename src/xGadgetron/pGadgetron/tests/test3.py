'''Test set 3.

Usage:
  test3 [--help | options]

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

import math
import sys

from pGadgetron import *

def test_main(rec = False, verb = False, throw = True):

    test = pTest('test3.txt', rec, throw = throw)
    test.verbose = verb

    data_path = mr_data_path()
    input_data = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')
    input_norm = input_data.norm()
    test.check(input_norm)
    alt_norm = math.sqrt(abs(input_data*input_data))
    test.check(abs(alt_norm/input_norm - 1), abs_tol = 1e-4)

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    processed_norm = processed_data.norm()
    test.check(processed_norm)

    for i in range(2):
        acq = processed_data.acquisition(i)
        print('--- acquisition %d' % i)
        for p in [ \
            'flags', 'kspace_encode_step_1', \
            'slice', 'repetition']:
            form = p + ' %d'
            test.check(acq.info(p))
            #print(form % acq.info(p))

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()
    images_norm = complex_images.norm()
    test.check(images_norm)
    alt_norm = math.sqrt(abs(complex_images*complex_images))
    test.check(abs(alt_norm/images_norm - 1), abs_tol = 1e-4)

    for i in range(complex_images.number()):
        complex_image = complex_images.image(i)
        print('--- image %d' % i)
        for p in [ \
            'version', 'flags', 'data_type', 'channels', \
            'slice', 'repetition', \
            'image_type', 'image_index', 'image_series_index' \
            ]:
            form = p + ' %d'
            test.check(complex_image.info(p))
            #print(form % complex_image.info(p))

    return test.failed, test.ntest

if __name__ == '__main__':

    __version__ = '0.1.0'
    from docopt import docopt
    args = docopt(__doc__, version=__version__)

    record = args['--record']
    verbose = args['--verbose']

    try:
        failed, ntest = test_main(record, verbose, throw = False)
        if failed == 0:
            if not record:
                print('all %d tests passed' % ntest)
            else:
                print('%d measurements recorded' % ntest)
            sys.exit(0)
        else:
            print('%d of %d tests failed' % (failed, ntest))
            sys.exit(failed)

    except error as err:
        # display error information
        print('??? %s' % err.value)
        sys.exit(-1)
