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
import os
from sirf.Utilities import runner, pTest, RE_PYEXT, __license__
from sirf.Utilities import examples_data_path
from sirf.Utilities import data_container_algebra_tests
#from sirf.Utilities import test_data_container_algebra
from sirf.Reg import MessageRedirector
import sirf.Reg as reg

__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov"


def test_main(rec=False, verb=False, throw=True, no_ret_val=True):
    MessageRedirector()

    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb
    data_path = examples_data_path('Registration')
    image = reg.ImageData(os.path.join(data_path, 'test2.nii.gz'))
    image /= image.norm()
    data_container_algebra_tests(test, image)
    #test_data_container_algebra(test, image)

    #return test.failed, test.ntest
    if no_ret_val:
        return
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__, no_ret_val=False)
