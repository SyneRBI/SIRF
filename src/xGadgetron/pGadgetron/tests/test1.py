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
import numpy
__version__ = "3.1.0"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    data_path = examples_data_path('MR')
    AcquisitionData.set_storage_scheme('memory')
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
    am = AcquisitionModel(processed_data, csms)
    am.set_coil_sensitivity_maps(csms)
    fwd_acqs = am.forward(complex_images)
    fwd_acqs_norm = fwd_acqs.norm()
    test.check(fwd_acqs_norm)

    rng = am.range_geometry()
    dom = am.domain_geometry()
    rng = rng - processed_data
    dom = dom - csms
    test.check_if_equal(0, rng.norm())
    test.check_if_equal(0, dom.norm())

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

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
