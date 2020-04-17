# -*- coding: utf-8 -*-
"""sirf.STIR tests
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
from sirf.STIR import *
from sirf.Utilities import runner, RE_PYEXT, __license__
__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


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

    # create an acq_model that is explicitly a RayTracingMatrix and test it a tiny bit
    am = AcquisitionModelUsingRayTracingMatrix()
    am.set_num_tangential_LORs(3);
    test.check_if_equal(3, am.get_num_tangential_LORs());

    # create the matrix on its own, and use that for later tests
    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)
    test.check_if_equal(2, matrix.get_num_tangential_LORs());

    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    data_path = examples_data_path('PET')

    raw_data_file = existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
    ad = AcquisitionData(raw_data_file)
    adata = ad.as_array()
    s = norm(adata)
    v = var(adata)
    test.check(s)
    test.check(v)

    filter = TruncateToCylinderProcessor()

    image_size = (31, 111, 111)
    voxel_size = (3.375, 3, 3)
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    test.check_if_equal(voxel_size, image.voxel_sizes())

    filter.apply(image)
    image_arr = image.as_array()
    s = norm(image_arr)
    v = var(image_arr)
    test.check(s)
    test.check(v)

    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.5)
    prior.set_up(image)

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

    # Test geom info
    geom_info = image.get_geometrical_info()
    geom_info.print_info()
    if geom_info.get_size() != (image_size[2],image_size[1],image_size[0]):
        raise AssertionError("SIRF get_geometrical_info().get_size() failed.")
    if geom_info.get_spacing() != (voxel_size[2],voxel_size[1],voxel_size[0]):
        raise AssertionError("SIRF get_geometrical_info().get_spacing() failed.")

    # Test zoom_image
    new_size = (3,2,5)
    zoomed_im = image.zoom_image(size=new_size)
    if zoomed_im.dimensions() != new_size:
        raise AssertionError("STIRImageData zoom_image() failed.\n\t" + \
            "Expected new size: " + str(new_size) + "\n\t" + \
            "Actual new size: " + str(zoomed_im.dimensions()))

    # Test move to scanner centre
    moved_im = image.move_to_scanner_centre(ad)

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
