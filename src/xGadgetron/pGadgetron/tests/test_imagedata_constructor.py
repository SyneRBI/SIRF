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
    print("Running the constructor test")

    data_path = examples_data_path('MR')
    AcquisitionData.set_storage_scheme('memory')
    
    rawdata = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')
    
    imgdata = ImageData()
    imgdata.from_acquisition_data(rawdata)

    test_successful = (imgdata.number() == rawdata.number())

    return test_successful, 1


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
