# -*- coding: utf-8 -*-
"""sirf.Gadgetron Test set 1.
v{version}

Fully sampled data tests

Usage:
  test1 [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
# Created on Tue Nov 21 10:17:28 2017
from sirf.Gadgetron import *
from sirf.Utilities import runner, RE_PYEXT, __license__
__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    data_path = examples_data_path('MR')
    input_data = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')
    test.check(input_data.norm())

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    test.check(processed_data.norm())

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()
    test.check(complex_images.norm())

    csms = CoilSensitivityData()

    processed_data.sort()
    csms.calculate(processed_data)
    am = AcquisitionModel(processed_data, complex_images)
    am.set_coil_sensitivity_maps(csms)
    fwd_acqs = am.forward(complex_images)
    fwd_acqs_norm = fwd_acqs.norm()
    test.check(fwd_acqs_norm)

    acqs_diff = fwd_acqs - processed_data
    rr = acqs_diff.norm()/fwd_acqs_norm
    test.check(rr, abs_tol = 1e-4)

    bwd_images = am.backward(processed_data)
    imgs_diff = bwd_images - complex_images
    rd = imgs_diff.norm()/complex_images.norm()
    test.check(rd, abs_tol = 1e-4)

##    xFy = processed_data * fwd_acqs
##    Bxy = bwd_images * complex_images
    xFy = processed_data.dot(fwd_acqs)
    Bxy = bwd_images.dot(complex_images)
    test.check(abs(xFy.real/Bxy.real - 1), abs_tol = 1e-4)
    test.check(abs(xFy.imag/xFy.real), abs_tol = 1e-4)
    test.check(abs(Bxy.imag/Bxy.real), abs_tol = 1e-4)

    pad2 = processed_data - processed_data
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr)
    print('acquisitions subtraction error: %.1e' % d)
    test.check_if_equal(0, d)
    processed_data.subtract(processed_data, out=pad2)
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr)
    print('acquisitions subtraction (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    pad2 = processed_data * processed_data
    pad_arr = processed_data.as_array()
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr*pad_arr)
    print('acquisitions multiplication error: %.1e' % d)
    test.check_if_equal(0, d)
    processed_data.multiply(processed_data, out=pad2)
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr*pad_arr)
    print('acquisitions multiplication (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    acq = processed_data.copy()
    pad2_arr[:] = 2.0
    acq.fill(pad2_arr)
    pad2 = processed_data / acq
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr/2)
    print('acquisitions division error: %.1e' % d)
    test.check_if_equal(0, d)
    processed_data.divide(acq, out=pad2)
    pad2_arr = pad2.as_array()
    d = numpy.linalg.norm(pad2_arr - pad_arr/2)
    print('acquisitions division (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    # testing in-place algebra
    acq = processed_data
    acq_clone = acq.clone()
    acq_clone /= acq_clone
    acq_clone *= acq
    acq_clone -= acq
    d = acq_clone.norm()/acq.norm()
    print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)
    acq_clone += acq
    d = (acq_clone.norm() - acq.norm())/acq.norm()
    print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)

    ci2 = complex_images - complex_images
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr)
    print('images subtraction error: %.1e' % d)
    test.check_if_equal(0, d)
    complex_images.subtract(complex_images, out=ci2)
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr)
    print('images subtraction (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    ci2 = complex_images * complex_images
    ci_arr = complex_images.as_array()
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr*ci_arr)
    print('images multiplication error: %.1e' % d)
    test.check_if_equal(0, d)
    complex_images.multiply(complex_images, out=ci2)
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr*ci_arr)
    print('images multiplication (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    img = complex_images.copy()
    ci2_arr[:] = 2.0
    img.fill(ci2_arr)
    ci2 = complex_images / img
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr/2)
    print('images division error: %.1e' % d)
    test.check_if_equal(0, d)
    complex_images.divide(img, out=ci2)
    ci2_arr = ci2.as_array()
    d = numpy.linalg.norm(ci2_arr - ci_arr/2)
    print('images division (with out=) error: %.1e' % d)
    test.check_if_equal(0, d)
    images_copy = complex_images.copy()
    # testing /=, *= and -=
    images_copy /= images_copy
    images_copy *= complex_images
    images_copy -= complex_images
    d = images_copy.norm()/complex_images.norm()
    print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)
    images_copy += complex_images
    d = (images_copy.norm() - complex_images.norm())/complex_images.norm()
    print('%f is 0.0' % d)
    test.check_if_equal(1, d < 1e-6)

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
