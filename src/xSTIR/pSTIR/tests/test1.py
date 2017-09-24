''' pSTIR tests
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

from pSTIR import *

def test_failed(ntest, expected, actual, abstol, reltol):
    if abs(expected - actual) < abstol + reltol*expected:
        print('+++ test %d passed' % ntest)
        return 0
    else:
        print('+++ test %d failed' % ntest)
        return 1

def norm(v):
    vv = v*v
    nv = v.size
#    return vv.sum()/nv
    return math.sqrt(vv.sum()/nv)

# a function to compute the variance after conversion to double to avoid
# rounding problems with older numpy versions
def var(v):
    return v.astype(numpy.float64).var()

def main():

    failed = 0
    eps = 1e-4

    # create matrix to be used by the acquisition model
    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # locate the input data file folder
    data_path = petmr_data_path('pet')

    # PET acquisition data to be read from this file
    raw_data_file = existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
    ad = AcquisitionData(raw_data_file)
    adata = ad.as_array()
    s = norm(adata)
    v = var(adata)
    failed += test_failed(1, 2.510818, s, 0, eps)
    failed += test_failed(2, 5.444323, v, 0, eps)
    #print('acquisitions mean sum of squares: %f, variance: %f' % (s, v))

    # create filter
    filter = TruncateToCylinderProcessor()

    # create initial image estimate
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375) # voxel sizes are in mm
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    
    filter.apply(image)
    image_arr = image.as_array()
    s = norm(image_arr)
    v = var(image_arr)
    failed += test_failed(3, 0.876471, s, 0, eps)
    failed += test_failed(4, 0.178068, v, 0, eps)
    #print('image mean sum of squares: %f, variance: %f' % (s, v))

    # create prior
    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.5)

    # set number of subsets
    num_subsets = 12

    # create objective function
    obj_fun = make_Poisson_loglikelihood(ad)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_num_subsets(num_subsets)
    obj_fun.set_up(image)

    # select subset
    subset = 0

    # get sensitivity as ImageData
    ss_img = obj_fun.get_subset_sensitivity(subset)

    # get back projection of the ratio of measured to estimated acquisition data
    grad_img = obj_fun.get_backprojection_of_acquisition_ratio(image, subset)

    # get gradient of prior as ImageData
    pgrad_img = prior.get_gradient(image)

    # copy to Python arrays
    image_arr = image.as_array()
    ss_arr = ss_img.as_array()
    grad_arr = grad_img.as_array()
    pgrad_arr = pgrad_img.as_array()

    # update image data
    ss_arr[ss_arr < 1e-6] = 1e-6 # avoid division by zero
    update = grad_arr/(ss_arr + pgrad_arr/num_subsets)
    image_arr = image_arr*update

    s = norm(image_arr)
    v = var(image_arr)
    failed += test_failed(5, 0.012314, s, 0, eps)
    failed += test_failed(6, 0.000052, v, eps, eps)
    #print('image mean sum of squares: %f, variance: %f' % (s, v))
    s = norm(update)
    v = var(update)
    failed += test_failed(7, 3.846513, s, 0, eps)
    failed += test_failed(8, 14.775219, v, 0, eps)
    #print('update mean sum of squares: %f, variance: %f' % (s, v))
    s = norm(ss_arr)
    v = var(ss_arr)
    failed += test_failed(9, 27.990159, s, 0, eps)
    failed += test_failed(10, 207.401144, v, 0, eps)
    #print('sensitivity mean sum of squares: %f, variance: %f' % (s, v))
    s = norm(grad_arr)
    v = var(grad_arr)
    failed += test_failed(11, 98.049032, s, 0, eps)
    failed += test_failed(12, 9599.796540, v, 0, eps)
    #print('gradient mean sum of squares: %f, variance: %f' % (s, v))
    s = norm(pgrad_arr)
    v = var(pgrad_arr)
    failed += test_failed(13, 0.710633, s, 0, eps)
    failed += test_failed(14, 0.505000, v, 0, eps)
    #print('prior gradient mean sum of squares: %f, variance: %f' % (s, v))
    return failed

try:
    failed = main()
    if failed == 0:
        print('all tests passed')
        sys.exit(0)
    else:
        print('%d tests failed' % failed)
        sys.exit(failed)

except error as err:
    # display error information
    print('??? %s' % err.value)
    sys.exit(-1)



