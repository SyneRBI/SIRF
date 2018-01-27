# -*- coding: utf-8 -*-
"""pSTIR Acquisitions and images algebra tests
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
from pSTIR import *
__version__ = "0.2.2"
__author__ = "Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    msg_red = MessageRedirector()

    data_path = petmr_data_path('pet')
    raw_data_file = existing_filepath(data_path, 'my_forward_projection.hs')
    acq_data = AcquisitionData(raw_data_file)

    if verb:
        print('Checking acquisition data algebra:')
    new_acq_data = acq_data.clone()
    diff = new_acq_data - acq_data
    test.check(diff.norm())
    test.check(1 - math.sqrt(acq_data * acq_data) / acq_data.norm())
    new_acq_data = acq_data * 10.0
    test.check(1 - 10 * acq_data.norm() / new_acq_data.norm())

    if verb:
        print('Checking images algebra:')
    image_data = acq_data.create_uniform_image(10.0)
    diff = image_data.clone() - image_data
    test.check(diff.norm())
    test.check(1 - math.sqrt(image_data * image_data) / image_data.norm())
    new_image_data = image_data * 10
    test.check(1 - 10 * image_data.norm() / new_image_data.norm())

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
