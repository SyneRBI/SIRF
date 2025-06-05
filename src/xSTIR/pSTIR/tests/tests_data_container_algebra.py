# -*- coding: utf-8 -*-
"""sirf.STIR algebra tests
v{version}

Usage:
  tests_data_container_algebra [--help | options]

Options:
  -r, --record   not used
  -v, --verbose  report each test status

{author}

{licence}
"""
import numpy
from sirf.STIR import *
from sirf.Utilities import runner, RE_PYEXT, __license__
from sirf.Utilities import data_container_algebra_tests

__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov"


def test_main(rec=False, verb=False, throw=True):
    MessageRedirector()

    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb
    data_path = examples_data_path('PET')

    for scheme in ("file", "memory"):
        AcquisitionData.set_storage_scheme(scheme)
        raw_data_file = existing_filepath(data_path, 'Utahscat600k_ca_seg4.hs')
        ad = AcquisitionData(raw_data_file)
        data_container_algebra_tests(test, ad)

        image_size = (31, 111, 111)
        voxel_size = (3.375, 3, 3)
        image = ImageData()
        image.initialise(image_size, voxel_size)
        image.fill(1.0)
        data_container_algebra_tests(test, image)

    #return test.failed, test.ntest
    numpy.testing.assert_equal(test.failed, 0)


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
