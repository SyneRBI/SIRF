# -*- coding: utf-8 -*-
"""sirf.Gadgetron test set 2.
v{version}

Undersampled data tests

Usage:
  test2 [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
# Created on Tue Nov 21 11:23:39 2017
from sirf.Gadgetron import *
from sirf.Utilities import is_operator_adjoint, runner, RE_PYEXT, __license__
__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    data_path = examples_data_path('MR')
    input_data = AcquisitionData(
            data_path + '/simulated_MR_2D_cartesian_Grappa2.h5')
    test.check(input_data.norm())

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    test.check(processed_data.norm())

    recon = CartesianGRAPPAReconstructor()
    recon.compute_gfactors(False)
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

    if not is_operator_adjoint(am, max_err = 1e-3):
      raise AssertionError("Gadgetron operator is not adjoint")

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
