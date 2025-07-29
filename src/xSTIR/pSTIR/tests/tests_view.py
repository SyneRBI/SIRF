# -*- coding: utf-8 -*-
"""sirf.STIR data containers views tests
v{version}

Usage:
  tests_three [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import numpy
from sirf.SIRF import norm, dot, copyto
from sirf.STIR import MessageRedirector, AcquisitionData, ImageData
from sirf.Utilities import pTest, runner, RE_PYEXT, examples_data_path, existing_filepath
__version__ = "0.2.4"
__author__ = "Evgueni Ovtchinnikov, Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True, no_ret_val=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    _ = MessageRedirector(warn=None)
    AcquisitionData.set_storage_scheme('memory')

    data_path = examples_data_path('PET')
    raw_data_file = existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
    acq_data = AcquisitionData(raw_data_file)
    x = acq_data + 0
    y = x + 0
    z = x + 0

    x_view = x.asarray(copy=False)
    y_view = y.asarray(copy=False)
    z_view = z.asarray(copy=False)

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

    init_image_file = existing_filepath(data_path, 'test_image_PM_QP_6.hv')
    image_data = ImageData(init_image_file)
    x = image_data + 0
    y = x + 0
    z = x + 0

    x_view = x.asarray(copy=False)
    y_view = y.asarray(copy=False)
    z_view = z.asarray(copy=False)

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
