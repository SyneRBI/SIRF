
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


def test_2d_slice_stack(rec=False, verb=False, throw=True):

    print("Running the ImageData from AcquisitionData constructor test")
    data_path = '/media/sf_CCPPETMR/TestData/Input/xGadgetron/pGadgetron'

    rawdata = AcquisitionData(data_path + '/CV_2D_Stack_144.h5')
    rawdata = preprocess_acquisition_data(rawdata)
    rawdata.sort()
    
    imgdata = ImageData()
    imgdata.from_acquisition_data(rawdata)
    
    print(imgdata.get_ISMRMRD_info('position'))

    print("The rawdata directions are: {}".format(rawdata.get_ISMRMRD_info('read_dir')[0]))
    print("The rawdata directions are: {}".format(rawdata.get_ISMRMRD_info('phase_dir')[0]))
    print("The rawdata directions are: {}".format(rawdata.get_ISMRMRD_info('slice_dir')[0]))

    print("The image data directions are: {}".format(imgdata.get_ISMRMRD_info('read_dir')[0]))
    print("The image data directions are: {}".format(imgdata.get_ISMRMRD_info('phase_dir')[0]))
    print("The image data directions are: {}".format(imgdata.get_ISMRMRD_info('slice_dir')[0]))

    # test_failed = not test_successful
    test_failed = True
    return test_failed, 1


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
