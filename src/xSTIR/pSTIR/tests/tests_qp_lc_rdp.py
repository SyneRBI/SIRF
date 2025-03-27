# -*- coding: utf-8 -*-
"""sirf.STIR prior tests (including stir.STIR.QuadraticPrior, stir.STIR.LogcoshPrior 
and sirf.STIR.RelativeDifferencePrior tests)
v{version}

Usage:
  tests_qp_lc_rdp [--help | options]

Options:
  -r, --record   record the measurements rather than check them
  -v, --verbose  report each test status

{author}

{licence}
"""
import sirf.STIR
import sirf.config
from sirf.Utilities import runner, RE_PYEXT, __license__, examples_data_path, pTest
__version__ = "2.1.0"
__author__ = "Imraj Singh, Evgueni Ovtchinnikov, Kris Thielemans, Edoardo Pasca"

try:
    subprocess.check_output('nvidia-smi')
    has_nvidia = True
except:
    has_nvidia = False


def Hessian_test(test, prior, x, eps=1e-3):
    """Checks that grad(x + dx) - grad(x) is close to H(x)*dx
        """
    if x.norm() > 0:
        dx = x.clone()
        dx *= eps/dx.norm()
        dx += eps/2
    else:
        dx = x + eps
    y = x + dx
    gx = prior.gradient(x.clone())
    gy = prior.gradient(y.clone())
    dg = gy - gx
    Hdx = prior.multiply_with_Hessian(x, dx)
    gynorm = gy.norm()
    if gynorm > 0:
        q = (dg - Hdx).norm()/gynorm
    else:
        q = (dg - Hdx).norm()
    #print('norm of grad(x + dx) - grad(x): %f' % dg.norm())
    #print('norm of H(x)*dx: %f' % Hdx.norm())
    #print('relative difference: %g' % q)
    test.check_if_less(q, .01*eps)


def test_main(rec=False, verb=False, throw=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    _ = sirf.STIR.MessageRedirector(warn=None)

    im_0 = sirf.STIR.ImageData(examples_data_path('PET') + '/brain/emission.hv')
    im_1 = im_0.get_uniform_copy(1)
    im_2 = im_0.get_uniform_copy(2)

    for im in [im_0, im_1, im_2]:
      for penalisation_factor in [0,1,4]:
        for kappa in [True, False]:
          priors = [sirf.STIR.QuadraticPrior(), sirf.STIR.LogcoshPrior(), sirf.STIR.RelativeDifferencePrior()]
          if sirf.config.STIR_WITH_CUDA and has_nvidia:
              priors.append(sirf.STIR.CudaRelativeDifferencePrior())
          for prior in priors:
            if kappa:
              prior.set_kappa(im_0)
              # Check if the kappa is stored/returned correctly
              test.check_if_equal_within_tolerance(im_0.norm(),prior.get_kappa().norm())

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
                   prior.set_kappa(im_0*2)
                   prior.set_up(im)
                   test.check_if_equal_within_tolerance(prior.get_gradient(im).norm(), grad_norm)

            if isinstance(prior, sirf.STIR.RelativeDifferencePrior):
                prior.set_epsilon(im.max()*.01)
            Hessian_test(test, prior, im, 0.03)
            
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__)
