'''pSTIR tests

Usage:
  tests_one [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
2017 Casper da Costa-Luis

This is software developed for the Collaborative Computational
Project in Positron Emission Tomography and Magnetic Resonance imaging
(http://www.ccppetmr.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at
      http://www.apache.org/licenses/LICENSE-2.0
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
'''
import math
from pSTIR import *
from os import path

def norm(v):
    vv = v*v
    nv = v.size
    return math.sqrt(vv.sum()/nv)

def var(v):
    """function to compute the variance after conversion to double to avoid
    rounding problems with older numpy versions
    """
    return v.astype(numpy.float64).var()

def test_main(rec=False, verb=False, throw=True, eps=1e-4):

    msg_red = MessageRedirector()

    datafile = path.join(path.dirname(__file__), 'test1.txt')
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    data_path = petmr_data_path('pet')

    raw_data_file = existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
    ad = AcquisitionData(raw_data_file)
    adata = ad.as_array()
    s = norm(adata)
    v = var(adata)
    test.check(s)
    test.check(v)

    filter = TruncateToCylinderProcessor()

    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375)
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    
    filter.apply(image)
    image_arr = image.as_array()
    s = norm(image_arr)
    v = var(image_arr)
    test.check(s)
    test.check(v)

    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.5)

    num_subsets = 12

    obj_fun = make_Poisson_loglikelihood(ad)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_num_subsets(num_subsets)
    if verb:
        print('setting up objective function, please wait...')
    obj_fun.set_up(image)

    subset = 0

    ss_img = obj_fun.get_subset_sensitivity(subset)

    grad_img = obj_fun.get_backprojection_of_acquisition_ratio(image, subset)

    pgrad_img = prior.get_gradient(image)

    image_arr = image.as_array()
    ss_arr = ss_img.as_array()
    grad_arr = grad_img.as_array()
    pgrad_arr = pgrad_img.as_array()

    ss_arr[ss_arr < 1e-6] = 1e-6
    update = grad_arr/(ss_arr + pgrad_arr/num_subsets)
    image_arr = image_arr*update

    s = norm(image_arr)
    v = var(image_arr)
    test.check(s)
    test.check(v)
    s = norm(update)
    v = var(update)
    test.check(s)
    test.check(v)
    s = norm(ss_arr)
    v = var(ss_arr)
    test.check(s)
    test.check(v)
    s = norm(grad_arr)
    v = var(grad_arr)
    test.check(s)
    test.check(v)
    s = norm(pgrad_arr)
    v = var(pgrad_arr)
    test.check(s)
    test.check(v)

    return test.failed, test.ntest

if __name__ == '__main__':

    __version__ = '0.2.0'
    from docopt import docopt
    args = docopt(__doc__, version=__version__)
    record = args['--record']
    verbose = args['--verbose']

    failed, ntest = test_main(record, verbose, throw=False)
    if failed:
        print('%d/%d tests failed' % (failed, ntest))
        sys.exit(failed)
    print('all %d tests passed' % ntest)
