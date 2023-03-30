# -*- coding: utf-8 -*-
"""sirf.Gadgetron algebra tests
v{version}

Usage:
  test5 [--help | options]

Options:
  -r, --record   not used
  -v, --verbose  report each test status

{author}

{licence}
"""
import math
import warnings
from sirf.Gadgetron import *
from sirf.Utilities import runner, RE_PYEXT, __license__
from sirf.Utilities import test_data_container_algebra

__version__ = "0.2.3"
__author__ = "Evgueni Ovtchinnikov"


def test_main(rec=False, verb=False, throw=True):

    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb
    data_path = examples_data_path('MR')

    input_data = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    test_data_container_algebra(test, processed_data)

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()
    test_data_container_algebra(test, complex_images)

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
