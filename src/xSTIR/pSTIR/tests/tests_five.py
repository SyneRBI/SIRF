# -*- coding: utf-8 -*-
"""STIR DataContainer algebra tests
v{version}

Usage:
  tests_five [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
from sirf.Utilities import petmr_data_path, existing_filepath, \
                           pTest, RE_PYEXT , runner
from sirf.STIR import MessageRedirector, AcquisitionData
__version__ = "0.2.2"
__author__ = "Casper da Costa-Luis, Edoardo Pasca"


def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = True

    msg_red = MessageRedirector()

    data_path = petmr_data_path('pet')
    raw_data_file = existing_filepath(data_path, 'my_forward_projection.hs')
    acq_data = AcquisitionData(raw_data_file)

    if verb:
        print('Checking images algebra:')
    image_data = acq_data.create_uniform_image(1.0)
    dims = image_data.dimensions()
    # N number of elements in the array
    N = 1
    for i,el in enumerate(dims):
        N *= el
    
    # 1 test sum: N * 1 / N = 1
    test.check(image_data.sum()/N)
    # test algebra 2 to 5
    # 2 DataContainer add (2+1) = 3
    # image_data = acq_data.create_uniform_image(1.0)
    b = acq_data.create_uniform_image(2.0)
    c = b + image_data
    test.check(c.sum()/N)
    # 3 DataContainer subtract 1 - (2) = -1
    image_data = acq_data.create_uniform_image(1.0)
    b = acq_data.create_uniform_image(2.0)
    c = image_data - b
    test.check(c.sum()/N)
    # 4 DataContainer multiply ( 2 * 1 ) = 2
    c = b * image_data
    test.check(c.sum()/N)
    # 5 DataContainer divide (1 / 2) = 0.5
    c = image_data.divide(b)
    test.check(c.sum()/N)
    # 6 power
    b = 1.5*image_data
    c = b.power(0.5)
    test.check(c.sum()/N)
    # 7 maximum
    test.check(c.maximum(b).sum()/N)
    # 8 sign
    b = -1 * image_data
    test.check(b.sign().sum()/N)
    # 9 abs
    test.check(b.abs().sum()/N)
    # 10 sqrt
    b = 1.5*image_data
    c = b.sqrt()
    test.check(c.sum()/N)
    # inline algebra
    # 11 inline add [1.5] + [1] = [2.5]
    b = acq_data.create_uniform_image(1.5)
    b += image_data
    test.check(b.sum()/N)
    # 12 inline subtract [1.5] - [1]
    b = acq_data.create_uniform_image(1.5)
    b -= image_data
    test.check(b.sum()/N)
    # 13 inline multiply
    b = acq_data.create_uniform_image(1.8)
    b *= 2.5
    # b = b.imul(2.5)
    test.check(b.sum()/N)
    # 14 inline divide
    b = acq_data.create_uniform_image(1.5)
    b /= 3
    test.check(b.sum()/N)

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
