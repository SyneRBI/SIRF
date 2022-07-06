# -*- coding: utf-8 -*-
"""sirf.STIR.QuadraticPrior and sirf.STIR.RelativeDifferencePrior tests
v{version}

Usage:
  tests_qp_rdp [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import math
import sirf.STIR as pet
from sirf.Utilities import runner, RE_PYEXT, __license__, examples_data_path, pTest
import numpy
__version__ = "3.4.0"
__author__ = "Imraj Singh"


def test_kappa(prior, image, test):
    prior.set_penalisation_factor(1)
    prior.set_kappa(image)
    prior.set_up(image)
    test.check_if_equal_within_tolerance(image.norm(),prior.get_kappa().norm())
    pgrad_img = prior.get_gradient(image)
    test.check(pgrad_img.norm())
    return test

def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    msg_red = pet.MessageRedirector()

    image = pet.ImageData(examples_data_path('PET') + '/thorax_single_slice/emission.hv')

    test = test_kappa(pet.RelativeDifferencePrior(), image, test)
    test = test_kappa(pet.QuadraticPrior(), image, test)

    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
