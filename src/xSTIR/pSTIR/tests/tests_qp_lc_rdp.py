# -*- coding: utf-8 -*-
"""sirf.STIR prior tests (including stir.STIR.QuadraticPrior, stir.STIR.LogcoshPrior 
and sirf.STIR.RelativeDifferencePrior tests)
v{version}

Usage:
  tests_qp_rdp [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
from sirf.STIR import *
from sirf.Utilities import runner, RE_PYEXT, __license__, examples_data_path, pTest
__version__ = "1.0.0"
__author__ = "Imraj Singh"
  

def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    msg_red = MessageRedirector(warn=None)

    im_thorax = ImageData(examples_data_path('PET') + '/thorax_single_slice/emission.hv')
    im_1 = im_thorax.get_uniform_copy(1)
    im_2 = im_thorax.get_uniform_copy(2)

    for im in [im_thorax, im_1, im_2]:
      for penalisation_factor in [0,1,4]:
        for kappa in [True, False]:
          for prior in [QuadraticPrior(), LogcoshPrior(), RelativeDifferencePrior()]:
            if kappa:
              prior.set_kappa(im_thorax)
              # Check if the kappa is stored/returned correctly
              test.check_if_equal_within_tolerance(im_thorax.norm(),prior.get_kappa().norm())

            prior.set_penalisation_factor(penalisation_factor)
            prior.set_up(im)
            # check if constant image or penalisation is zero
            if im.as_array().min()==im.as_array().max() or penalisation_factor == 0:
              # then grad norm should be zero
              test.check_if_equal_within_tolerance(prior.get_gradient(im).norm(),0)
            else:
              # or check against .txt
              grad_norm = prior.get_gradient(im).norm()
              test.check(grad_norm)
              if kappa:
                    # check if multiplying kappa and dividing penalisation factor gives same result
                   prior.set_penalisation_factor(penalisation_factor/4)
                   prior.set_kappa(im_thorax*2)
                   prior.set_up(im)
                   test.check_if_equal_within_tolerance(prior.get_gradient(im).norm(), grad_norm)
            
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
