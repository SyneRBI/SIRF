# -*- coding: utf-8 -*-
"""sirf.Gadgetron Views Test.
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
from sirf.SIRF import norm, dot, copyto
from sirf.Gadgetron import AcquisitionData, \
    AcquisitionDataView, ImageDataView, FullySampledReconstructor
from sirf.Utilities import pTest, runner, RE_PYEXT, examples_data_path
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
    x = processed_data
    y = x + 0
    z = x + 0

    x_view = AcquisitionDataView(x)
    y_view = AcquisitionDataView(y)
    z_view = AcquisitionDataView(z)

    y = 2*x
    norm_y = y.norm()
    copyto(y_view, x_view)
    y_view *= 2
    norm_y_view = norm(y_view)
    test.check_if_equal_within_tolerance(norm_y, norm_y_view, rel_tol=1e-5)
    z = x + y
    norm_z = z.norm()
    copyto(z_view, x_view)
    z_view += y_view
    norm_z_view = norm(z_view)
    test.check_if_equal_within_tolerance(norm_z, norm_z_view, rel_tol=1e-5)
    s = norm(x)
    t = norm(x_view)
    test.check_if_equal_within_tolerance(s, t, rel_tol=1e-5)
    s = x.sum()
    t = x_view.sum()
    test.check_if_equal_within_tolerance(s, t, rel_tol=1e-5)
    s = x.dot(y)
    t = dot(x_view, y_view)
    test.check_if_equal_within_tolerance(s, t, rel_tol=1e-5)

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    x = complex_images
    y = x + 0
    z = x + 0

    x_view = ImageDataView(x)
    y_view = ImageDataView(y)
    z_view = ImageDataView(z)

    y = 2*x
    norm_y = y.norm()
    copyto(y_view, x_view)
    y_view *= 2
    norm_y_view = norm(y_view)
    test.check_if_equal_within_tolerance(norm_y, norm_y_view, rel_tol=1e-5)
    z = x + y
    norm_z = z.norm()
    copyto(z_view, x_view)
    z_view += y_view
    norm_z_view = norm(z_view)
    test.check_if_equal_within_tolerance(norm_z, norm_z_view, rel_tol=1e-5)
    s = norm(x)
    t = norm(x_view)
    test.check_if_equal_within_tolerance(s, t, rel_tol=1e-5)
    s = x.sum()
    t = x_view.sum()
    test.check_if_equal_within_tolerance(s, t, rel_tol=1e-5)
    s = x.dot(y)
    t = dot(x_view, y_view)
    test.check_if_equal_within_tolerance(s, t, rel_tol=1e-5)

    numpy.testing.assert_equal(test.failed, 0)
    if no_ret_val:
        return
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__, no_ret_val=False)
