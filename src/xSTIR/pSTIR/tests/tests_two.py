# -*- coding: utf-8 -*-
"""sirf.STIR OSEM reconstruction tests
v{version}

Usage:
  tests_two [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import math
from sirf.STIR import *
from sirf.Utilities import runner, RE_PYEXT, __license__
import numpy
__version__ = "3.1.0"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"

def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    msg_red = MessageRedirector()

    for scheme in ("file", "memory"):
        AcquisitionData.set_storage_scheme(scheme)

        data_path = examples_data_path('PET')
        raw_data_file = existing_filepath(data_path, 'my_forward_projection.hs')
        acq_data = AcquisitionData(raw_data_file)
        test.check(acq_data.norm())

        image = acq_data.create_uniform_image(1.0)
        test.check_if_equal_within_tolerance(image.norm(), numpy.sqrt(numpy.prod(image.dimensions())))

        acq_model = AcquisitionModelUsingRayTracingMatrix()
        acq_model.set_up(acq_data, image)

        acq_mod_lin = acq_model.get_linear_acquisition_model()
        rng = acq_mod_lin.range_geometry()
        dom = acq_mod_lin.domain_geometry()
        rng = rng - acq_data
        dom = dom - image
        test.check_if_equal(0, rng.norm())
        test.check_if_equal(0, dom.norm())

        obj_fun = make_Poisson_loglikelihood(acq_data)
        obj_fun.set_acquisition_model(acq_model)

        num_subsets = 12
        recon = OSMAPOSLReconstructor()
        recon.set_objective_function(obj_fun)
        recon.set_num_subsets(num_subsets)
        recon.set_input(acq_data)
        recon.set_up(image)

        recon.set_current_estimate(image)

        num_iterations = 2
        for iteration in range(num_iterations):
            if verb:
                print('\n------------- iteration %d' % iteration)
            recon.update_current_estimate()
        test.check(recon.get_output().norm())

        if verb:
            print('projecting...')
        simulated_data = acq_model.forward(recon.get_output())
        diff = simulated_data * (
                acq_data.norm() / simulated_data.norm()) - acq_data
        res = diff.norm() / acq_data.norm()
        if verb:
            print('relative residual norm: %e' % res)
        test.check_if_zero_within_tolerance(res, abs_tol=0.28) # only 2 iterations, so low tolerance

        acq_copy = acq_data.get_uniform_copy(1.0)
        acq_copy *= acq_data
        diff = acq_copy
        diff -= acq_data
        d = diff.norm()
        print('in-place algebra error: %.1e' % d)
        test.check_if_equal(0, d)
        acq_copy += acq_data
        diff = acq_data - acq_copy
        diff_arr = diff.as_array()
        d = numpy.linalg.norm(diff_arr)
        print('acquisitions subtraction error: %.1e' % d)
        test.check_if_equal(0, d)
        acq_data.subtract(acq_data, out=diff)
        diff_arr = diff.as_array()
        d = numpy.linalg.norm(diff_arr)
        print('acquisitions subtraction (with out=) error: %.1e' % d)
        test.check_if_equal(0, d)
        ad2 = acq_data * acq_data
        ad_arr = acq_data.as_array()
        ad2_arr = ad2.as_array()
        d = numpy.linalg.norm(ad2_arr - ad_arr*ad_arr)
        print('acquisitions multiplication error: %.1e' % d)
        test.check_if_equal(0, d)
        acq_data.multiply(acq_data, out=ad2)
        ad2_arr = ad2.as_array()
        d = numpy.linalg.norm(ad2_arr - ad_arr*ad_arr)
        print('acquisitions multiplication (with out=) error: %.1e' % d)
        test.check_if_equal(0, d)
        ad2_arr[:] = 2.0
        acq_copy.fill(ad2_arr)
        ad2 = acq_data / acq_copy
        ad2_arr = ad2.as_array()
        d = numpy.linalg.norm(ad2_arr - ad_arr/2)
        print('acquisitions division error: %.1e' % d)
        test.check_if_equal(0, d)
        acq_data.divide(acq_copy, out=ad2)
        ad2_arr = ad2.as_array()
        d = numpy.linalg.norm(ad2_arr - ad_arr/2)
        print('acquisitions division (with out=) error: %.1e' % d)
        test.check_if_equal(0, d)
        acq_arr = acq_data.as_array()
        acq_conj = acq_data.conjugate()
        acq_arr2 = acq_data.as_array()
        d = numpy.linalg.norm(acq_arr - acq_arr2)
        test.check_if_equal(0, d)
        acq_arr_conj_sirf = acq_conj.as_array()
        acq_arr_conj_numpy = numpy.conjugate(acq_arr)
        d = numpy.linalg.norm(acq_arr_conj_sirf - acq_arr_conj_numpy)
        test.check_if_equal(0, d)
        acq_data.conjugate(out=acq_data)
        acq_arr_conj_sirf = acq_data.as_array()
        d = numpy.linalg.norm(acq_arr_conj_sirf - acq_arr_conj_numpy)
        test.check_if_equal(0, d)

        image_copy = image.get_uniform_copy(1.0)
        image_copy *= image
        diff = image_copy
        diff -= image
        d = diff.norm()
        print('in-place algebra error: %.1e' % d)
        test.check_if_equal(0, d)
        image_copy += image
        im2 = image - image_copy
        im2_arr = im2.as_array()
        d = numpy.linalg.norm(im2_arr)
        print('images subtraction error: %.1e' % d)
        test.check_if_equal(0, d)
        image.subtract(image_copy, out=im2)
        im2_arr = im2.as_array()
        d = numpy.linalg.norm(im2_arr)
        print('images subtraction (with out=) error: %.1e' % d)
        test.check_if_equal(0, d)
        im2 = image * image
        im_arr = image.as_array()
        im2_arr = im2.as_array()
        d = numpy.linalg.norm(im2_arr - im_arr * im_arr)
        print('images multiplication error: %.1e' % d)
        test.check_if_equal(0, d)
        image.multiply(image, out=im2)
        im2_arr = im2.as_array()
        d = numpy.linalg.norm(im2_arr - im_arr * im_arr)
        print('images multiplication (with out=) error: %.1e' % d)
        test.check_if_equal(0, d)
        im2_arr[:] = 2.0
        image_copy.fill(im2_arr)
        im2 = image / image_copy
        im2_arr = im2.as_array()
        d = numpy.linalg.norm(im2_arr - im_arr/2)
        print('images division error: %.1e' % d)
        test.check_if_equal(0, d)
        image.divide(image_copy, out=im2)
        im2_arr = im2.as_array()
        d = numpy.linalg.norm(im2_arr - im_arr/2)
        print('images division (with out=) error: %.1e' % d)
        test.check_if_equal(0, d)
        # test on fill with scalar
        types = (2,2.0,numpy.int32(2),numpy.int64(2))
        try:
            twos = types + (numpy.float128(2),)
        except AttributeError:
            twos = types
        for n in twos:
            image.fill(n)
            test.check_if_zero_within_tolerance((image-2).norm())
        img = image
        img_arr = img.as_array()
        img_conj = img.conjugate()
        img_arr2 = img.as_array()
        d = numpy.linalg.norm(img_arr - img_arr2)
        test.check_if_equal(0, d)
        img_arr_conj_sirf = img_conj.as_array()
        img_arr_conj_numpy = numpy.conjugate(img_arr)
        d = numpy.linalg.norm(img_arr_conj_sirf - img_arr_conj_numpy)
        test.check_if_equal(0, d)
        img.conjugate(out=img)
        img_arr_conj_sirf = img.as_array()
        d = numpy.linalg.norm(img_arr_conj_sirf - img_arr_conj_numpy)
        test.check_if_equal(0, d)

    #return test.failed, test.ntest
    numpy.testing.assert_equal(test.failed, 0)


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
