import argparse
import numpy
import pylab
import os
import sys
sys.path.append(os.environ.get('CSTIR_SRC') + '/../pSTIR')
import scipy
from scipy import optimize
import time

from pStir import *

parser = argparse.ArgumentParser(description = \
'''
Steepest descent demo
''')
parser.add_argument('-t', '--tau', action = 'store', type = float, \
    default = -1, \
    help = 'steepest descent step parameter, use a negative value to opt \
for the optimal value, default -1')
parser.add_argument('-s', '--steps', action = 'store', type = int, default = 3, \
    help = 'number of steepest descent steps, default 3')
parser.add_argument('-v', '--verbose', action = 'store_true', default = False, \
    help = 'if present and optimal steepest descent step opted for, \
prints optimal step search info') 
args = parser.parse_args()

def main():

    # direct all information printing to a file
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    # direct all warning printing to a file
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)
    # direct all error printing to stdout
    error_printer = printerTo('stdout', ERROR_CHANNEL)

    # create matrix to be used by the acquisition model
    matrix = RayTracingMatrix()
    matrix.set_num_tangential_LORs(2)

    # create acquisition model
    am = AcquisitionModelUsingMatrix()
    am.set_matrix(matrix)

    # define acquisition data
    ad = AcquisitionData('my_forward_projection.hs')
    # plot acquisition data
    adata = ad.as_array()
    print(adata.shape)
    pylab.figure(1)
    pylab.title('acquisitions')
    pylab.imshow(adata[10,:,:])
    print('close window to continue')
    pylab.colorbar()
    pylab.show()

    # create prior
    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.001)

    # create filter
    filter = CylindricFilter()

    # create initial image estimate
    nx = 111
    ny = 111
    nz = 31
    image_size = (nx, ny, nz)
    voxel_size = (3, 3, 3.375)
    image = Image()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    filter.apply(image)

    # create objective function
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_zero_seg0_end_planes(True)
    obj_fun.set_max_segment_num_to_process(3)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    obj_fun.set_prior(prior)

    # create OSMAPOSL reconstructor
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(12)
    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # plot the initial image
    idata = image.as_array()
    pylab.figure(1)
    pylab.title('initial image')
    pylab.imshow(idata[20,:,:])
    print('close window to continue')
    pylab.colorbar()
    pylab.show()

    print('computing initial objective function value...')
    print('objective function value: %e' % (obj_fun.value(image)))

    tau = args.tau
    steps = args.steps
    if args.verbose:
        disp = 3
        print('NOTE: below f(x) is the negative of the objective function value')
    else:
        disp = 0
    eps = 1e-6

    for iter in range(steps):

        # obtain gradient
        grad = obj_fun.gradient(image, 0)
        filter.apply(grad)
        gdata = grad.as_array()

        # find maximal steepest descent step parameter t in image + t*grad 
        # such that the new image remains positive
        max_image = idata.max()
        max_grad = abs(gdata).max()
        gdata[abs(gdata) < eps] = eps
        r = idata/gdata
        gdata[idata <= 0] = 0
        d = r[r < 0]
        if d.shape[0] > 0:
            maxstep = -d.max()
        else:
            maxstep = abs(tau)*max_image/max_grad

        if tau < 0:
            # find the optimal t
            fun = lambda t: -obj_fun.value(image.fill(idata + x*gdata))
            t = scipy.optimize.fminbound \
                (fun, 0, maxstep, xtol = 1e-4, maxfun = 3, disp = disp)
        else:
            # t is such that the relative change in image is not greater than tau
            t = tau*max_image/max_grad
            if t > maxstep:
                t = maxstep

        # perform steepest descent step
        print('step %d, max change in image %e' % (iter, t*max_grad))
        idata = idata + t*gdata

        # filter the new image
        image.fill(idata)
        filter.apply(image)
        idata = image.as_array()

        # plot the new image
        pylab.figure(iter + 1)
        pylab.title('image')
        pylab.imshow(idata[20,:,:])
        print('close window to continue...')
        pylab.colorbar()
        pylab.show()

        # quit if the new image has negative values
        min_image = idata.min()
        if min_image < -eps:
            print('image minimum is negative: %e' % min_image)
            break

    if tau > 0 or disp == 0:
        print('computing final objective function value...')
        print('objective function value: %e' % (obj_fun.value(image)))

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
