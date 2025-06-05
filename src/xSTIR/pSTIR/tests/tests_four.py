# -*- coding: utf-8 -*-
"""sirf.STIR Acquisitions and images algebra tests
v{version}

Usage:
  tests_four [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import math
import numpy
from sirf.STIR import *
from sirf.Utilities import runner, RE_PYEXT, __license__
__version__ = "0.2.4"
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

        if verb:
            print('Checking acquisition data algebra:')
        new_acq_data = acq_data.clone()
        diff = new_acq_data - acq_data
        acq_data_norm = acq_data.norm()
        test.check_if_zero_within_tolerance(diff.norm() / acq_data_norm)
        test.check_if_equal_within_tolerance(acq_data_norm, math.sqrt(acq_data.dot(acq_data)))
    #    test.check(1 - math.sqrt(acq_data * acq_data) / acq_data.norm())
        new_acq_data = acq_data * 10.0
        test.check_if_equal_within_tolerance(10 * acq_data_norm, new_acq_data.norm())

        if verb:
            print('Checking images algebra:')
        image_data = acq_data.create_uniform_image(10.0)
        diff = image_data.clone() - image_data
        image_data_norm = image_data.norm()
        test.check_if_zero_within_tolerance(diff.norm()/image_data_norm)
        test.check_if_equal_within_tolerance(image_data_norm, math.sqrt(image_data.dot(image_data)))
    #    test.check(1 - math.sqrt(image_data * image_data) / image_data.norm())
        new_image_data = image_data * 10
        test.check_if_equal_within_tolerance(10 * image_data_norm, new_image_data.norm())

    #return test.failed, test.ntest
    numpy.testing.assert_equal(test.failed, 0)


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
