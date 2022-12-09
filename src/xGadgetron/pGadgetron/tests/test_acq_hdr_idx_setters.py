# -*- coding: utf-8 -*-
"""sirf.Gadgetron test.
v{version}

Constructor of MR image data from MR acquisition data

Usage:
  test_imagedata_constructor [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status
{author}

{licence}
"""

from sirf.Gadgetron import *
from sirf.Utilities import runner, RE_PYEXT, __license__
__version__ = "0.2.3"
__author__ = "Johannes Mayer"


def test_main(rec=False, verb=False, throw=True):

    print("Running the ImageData from AcquisitionData constructor test")
    data_path = examples_data_path('MR')

    rawdata = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')

    encoding_limit = (10,100,50)
    encoding_name = "repetition" # any other string will cause failure since the encoding was not performed.
    rawdata.set_encoding_limit(encoding_name, encoding_limit)

    print(rawdata.get_header())

    test_successful = True

    test_failed = not test_successful
    return test_failed, 1


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
