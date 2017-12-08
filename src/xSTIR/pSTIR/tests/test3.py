'''pSTIR test set 3.

OSSPS reconstruction tests.

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

from pSTIR import *

def test_main(rec = False, verb = False, throw = True):

    test = pTest('test3.txt', rec, throw)
    test.verbose = verb

    msg_red = MessageRedirector()

    data_path = petmr_data_path('pet')
    raw_data_file = existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
    acq_data = AcquisitionData(raw_data_file)
    test.check(acq_data.norm())

    init_image_file = existing_filepath(data_path, 'test_image_PM_QP_6.hv')
    image_data = ImageData(init_image_file)
    test.check(image_data.norm())

    acq_model = AcquisitionModelUsingRayTracingMatrix()

    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_prior(QuadraticPrior().set_penalisation_factor(0.5))

    recon = OSSPSReconstructor()
    recon.set_num_subsets(4)
    recon.set_num_subiterations(2)
    recon.set_objective_function(obj_fun)
    recon.set_input(acq_data)
    print('setting up, please wait...')
    recon.set_up(image_data)
    recon.set_current_estimate(image_data)
    print('reconstructing, please wait...')
    recon.process()
    image_data = recon.get_output()
    test.check(image_data.norm())

    return test.failed, test.ntest

if __name__ == '__main__':

    __version__ = '0.1.0'
    from docopt import docopt
    args = docopt(__doc__, version = __version__)
    record = args['--record']
    verbose = args['--verbose']

    try:
        failed, ntest = test_main(record, verbose, throw = False)
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
