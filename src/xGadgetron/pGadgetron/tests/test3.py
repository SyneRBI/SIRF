# -*- coding: utf-8 -*-
"""Test set 3.
v{version}

Fully sampled data tests

Usage:
  test3 [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
from pGadgetron import *
import math
# Created on Tue Nov 21 10:17:28 2017
__version__ = "0.2.0"
__author__ = "Casper da Costa-Luis"


def test_main(rec=False, verb=False, throw=True):
    datafile = __file__.replace(".py", ".txt")
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    data_path = mr_data_path()
    input_data = AcquisitionData(data_path + '/simulated_MR_2D_cartesian.h5')
    input_norm = input_data.norm()
    test.check(input_norm)
    alt_norm = math.sqrt(abs(input_data*input_data))
    test.check(abs(alt_norm/input_norm - 1), abs_tol = 1e-4)

    prep_gadgets = ['RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)
    processed_norm = processed_data.norm()
    test.check(processed_norm)

    for i in range(2):
        acq = processed_data.acquisition(i)
        print('--- acquisition %d' % i)
        for p in [ \
            'flags', 'kspace_encode_step_1', \
            'slice', 'repetition']:
            form = p + ' %d'
            test.check(acq.info(p))
            #print(form % acq.info(p))

    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()
    images_norm = complex_images.norm()
    test.check(images_norm)
    alt_norm = math.sqrt(abs(complex_images*complex_images))
    test.check(abs(alt_norm/images_norm - 1), abs_tol = 1e-4)

    for i in range(complex_images.number()):
        complex_image = complex_images.image(i)
        print('--- image %d' % i)
        for p in [ \
            'version', 'flags', 'data_type', 'channels', \
            'slice', 'repetition', \
            'image_type', 'image_index', 'image_series_index' \
            ]:
            form = p + ' %d'
            test.check(complex_image.info(p))
            #print(form % complex_image.info(p))

    return test.failed, test.ntest


if __name__ == "__main__":
    test_runner(test_main, __doc__, __version__, __author__)
