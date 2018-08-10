# -*- coding: utf-8 -*-
"""pSTIR tests
v{version}

Usage:
  tests_one [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import math
from pSTIR import *
__version__ = "0.2.2"
__author__ = "Casper da Costa-Luis"


def norm(v):
    vv = v * v
    nv = v.size
    return math.sqrt(vv.sum() / nv)


def var(v):
    """function to compute the variance after conversion to double to avoid
    rounding problems with older numpy versions
    """
    from numpy import float64
    return v.astype(float64).var()


def test_main(rec=False, verb=False, throw=True):
    msg_red = MessageRedirector()

    datafile = RE_PYEXT.sub(".txt", __file__)
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
    prior.set_up(image)

    prior2 = PLSPrior()
    prior2.set_penalisation_factor(0.5)
    prior2.set_anatomical_image(image)
    prior2.set_up(image)

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
    update = grad_arr / (ss_arr + pgrad_arr / num_subsets)
    image_arr = image_arr * update

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


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
