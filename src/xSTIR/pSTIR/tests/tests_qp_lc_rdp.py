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
import subprocess
import numpy
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
        dx *= eps
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
    # print('norm of x: %f, dx: %f' % (x.norm(), dx.norm()))
    # print('norm of grad(x): %f, grad(x + dx): %f' % (gx.norm(), gy.norm()))
    # print('norm of grad(x + dx) - grad(x): %f' % dg.norm())
    # print('norm of H(x)*dx: %f' % Hdx.norm())
    # print('relative difference: %g' % q)
    if issubclass(type(prior), sirf.STIR.RelativeDifferencePrior):
       if x.min() == x.max() and dx.min() == dx.max():
          # skip test in this case, as grad(x+dx) = grad(x) = 0, but H dx is not (even analytically),
          # although it is small.
          # The difficult is knowing what the tolerance is for this test in that case
          q = 0
    test.check_if_less(q, .01*eps)


def Hessian_diagonal_test(test, prior, x, eps=1e-4):
    """Checks compute_Hessian_diagonal against numerical finite differences.

    For a subset of interior voxels i, checks:
        H_diag[i] ≈ (grad_i(x + eps*e_i) - grad_i(x)) / eps
    """
    H_diag = prior.compute_Hessian_diagonal(x)
    gx = prior.gradient(x.clone())

    x_arr = x.as_array()
    gx_arr = gx.as_array()
    H_diag_arr = H_diag.as_array()

    # Test a subset of interior voxels (avoid boundary where fewer neighbours)
    shape = x_arr.shape
    z_mid, y_mid, x_mid = shape[0] // 2, shape[1] // 2, shape[2] // 2
    test_coords = [(z_mid + dz, y_mid + dy, x_mid + dx)
                   for dz in [-1, 0, 1] for dy in [-1, 0, 1] for dx in [0]]

    for (iz, iy, ix) in test_coords:
        x_pert = x.clone()
        x_pert_arr = x_pert.as_array()
        x_pert_arr[iz, iy, ix] += eps
        x_pert.fill(x_pert_arr)

        g_pert = prior.gradient(x_pert)
        g_pert_arr = g_pert.as_array()

        numerical = (g_pert_arr[iz, iy, ix] - gx_arr[iz, iy, ix]) / eps
        analytical = H_diag_arr[iz, iy, ix]

        max_H = max(abs(analytical), abs(numerical), 1e-10)
        test.check_if_less(abs(numerical - analytical) / max_H, 0.01)


def test_main(rec=False, verb=False, throw=True, no_ret_val=True):
    datafile = RE_PYEXT.sub(".txt", __file__)
    test = pTest(datafile, rec, throw=throw)
    test.verbose = verb

    _ = sirf.STIR.MessageRedirector(warn=None)

    im_0 = sirf.STIR.ImageData(examples_data_path('PET') + '/brain/emission.hv')
    im_1 = im_0.get_uniform_copy(1)
    im_2 = im_0.get_uniform_copy(2)

    for im in [im_0, im_1, im_2]:
      # print('-------------- new image (see test source) -----------------')
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
            # print(f"penalisation_factor {penalisation_factor}, kappa {kappa}, prior {type(prior).__name__}")
            Hessian_test(test, prior, im, 0.03)

            if isinstance(prior, sirf.STIR.RelativeDifferencePrior) and penalisation_factor > 0:
                Hessian_diagonal_test(test, prior, im)

    numpy.testing.assert_equal(test.failed, 0)
    if no_ret_val:
        return
    return test.failed, test.ntest


if __name__ == "__main__":
    runner(test_main, __doc__, __version__, __author__, no_ret_val=False)
