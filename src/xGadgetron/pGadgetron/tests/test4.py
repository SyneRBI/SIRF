# -*- coding: utf-8 -*-
"""sirf.Gadgetron Test set 1.
v{version}

DataContainer algebra tests

Usage:
  test4 [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
from sirf.Gadgetron import *
from sirf.Utilities import runner, RE_PYEXT, __license__
import numpy
__version__ = "3.1.0"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True, no_ret_val=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    data_path = examples_data_path('MR')
    AcquisitionData.set_storage_scheme('memory')
    input_data = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    pad2 = processed_data - processed_data
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr)
    #print('acquisitions subtraction error: %.1e' % d)
    test.check_if_equal(0, d)
    processed_data.subtract(processed_data, out=pad2)
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr)
    #print('acquisitions subtraction (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    pad2 = processed_data * processed_data
    pad_arr = processed_data.as_array()
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr*pad_arr)
    #print('acquisitions multiplication error: %.1e' % d)
    test.check_if_equal(0, d)
    processed_data.multiply(processed_data, out=pad2)
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr*pad_arr)
    #print('acquisitions multiplication (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    acq = processed_data.copy()
    pad2_arr[:] = 2.0
    acq.fill(pad2_arr)
    pad2 = processed_data / acq
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr/2)
    #print('acquisitions division error: %.1e' % d)
    test.check_if_equal(0, d)
    processed_data.divide(acq, out=pad2)
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr/2)
    #print('acquisitions division (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    # testing in-place algebra
    acq = processed_data
    acq_clone = acq.clone()
    acq_clone /= acq_clone
    acq_clone *= acq
    acq_clone -= acq
    d = acq_clone.norm()/acq.norm()
    #print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)
    acq_clone += acq
    d = (acq_clone.norm() - acq.norm())/acq.norm()
    #print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)
    acq_arr = acq.as_array()
    acq_conj = acq.conjugate()
    acq_arr2 = acq.as_array()
    d = numpy.linalg.norm(acq_arr - acq_arr2)
    test.check_if_equal(0, d)
    acq_arr_conj_sirf = acq_conj.as_array()
    acq_arr_conj_numpy = numpy.conjugate(acq_arr)
    d = numpy.linalg.norm(acq_arr_conj_sirf - acq_arr_conj_numpy)
    test.check_if_equal(0, d)
    acq.conjugate(out=acq)
    acq_arr_conj_sirf = acq.as_array()
    d = numpy.linalg.norm(acq_arr_conj_sirf - acq_arr_conj_numpy)
    test.check_if_equal(0, d)

    ci_abs1 = numpy.abs(complex_images.as_array())
    ci_abs2 = complex_images.abs().as_array()
    d = numpy.linalg.norm(ci_abs1 - ci_abs2)
    test.check_if_equal(0, d)
    ci2 = ImageData()
    complex_images.abs(out=ci2)
    ci_abs3 = numpy.abs(ci2.as_array())
    d = numpy.linalg.norm(ci_abs1 - ci_abs3)
    test.check_if_equal(0, d)
    ci2 = complex_images - complex_images
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr)
    #print('images subtraction error: %.1e' % d)
    test.check_if_equal(0, d)
    complex_images.subtract(complex_images, out=ci2)
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr)
    #print('images subtraction (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    ci2 = complex_images * complex_images
    ci_arr = complex_images.as_array()
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr*ci_arr)
    #print('images multiplication error: %.1e' % d)
    test.check_if_equal(0, d)
    complex_images.multiply(complex_images, out=ci2)
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr*ci_arr)
    #print('images multiplication (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    img = complex_images.copy()
    ci2_arr[:] = 2.0
    img.fill(ci2_arr)
    ci2 = complex_images / img
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr/2)
    #print('images division error: %.1e' % d)
    test.check_if_equal(0, d)
    complex_images.divide(img, out=ci2)
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr/2)
    #print('images division (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    images_copy = complex_images.copy()
    # testing /=, *= and -=
    images_copy /= images_copy
    images_copy *= complex_images
    images_copy -= complex_images
    d = images_copy.norm()/complex_images.norm()
    #print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)
    images_copy += complex_images
    d = (images_copy.norm() - complex_images.norm())/complex_images.norm()
    #print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)
    # test on fill with scalar
    types = (2,2.0,numpy.int32(2),numpy.int64(2),complex(2))
    try:
        twos = types + (numpy.float128(2),)
    except AttributeError:
        twos = types
    for n in twos:
        complex_images.fill(n)
        test.check_if_zero_within_tolerance((complex_images-2).norm())
    img = complex_images
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

    numpy.testing.assert_equal(test.failed, 0)
    if no_ret_val:
        return
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__, no_ret_val=False)
