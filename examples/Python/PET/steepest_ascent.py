'''Steepest ascent demo.
Applies few steps of steepest ascent for the maximization of Poisson 
log-likelihood objective function using subset gradients.

Usage:
  steepest_ascent [--help | options]

Options:
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -s <nstp>, --steps=<nstp>   number of steepest descent steps [default: 3]
  -o, --optimal               use locally optimal steepest ascent
  -v, --verbose               verbose
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2020 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
## (http://www.ccpsynerbi.ac.uk/).
##
## Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##       http://www.apache.org/licenses/LICENSE-2.0
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')


# process command-line options
steps = int(args['--steps'])
opt = args['--optimal']
verbose = args['--verbose']
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
show_plot = not args['--non-interactive']


if opt:
    import scipy.optimize


def trunc(image):
    arr = image.as_array()
    arr[arr < 0 ] = 0
    out = image.copy()
    out.fill(arr)
    return out


def main():

    # engine's messages go to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # create acquisition model
    acq_model = AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # create filter that zeroes the image outside a cylinder of the same
    # diameter as the image xy-section size
    filter = TruncateToCylinderProcessor()

    # create initial image estimate
    nx = 111
    ny = 111
    nz = 31
    image_size = (nz, ny, nx)
    voxel_size = (3.375, 3, 3) # sizes are in mm
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    # apply the filter to the image
    filter.apply(image)

    # create objective function of Poisson logarithmic likelihood type
    # compatible with the acquisition data type
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_num_subsets(12)
    obj_fun.set_up(image)

    if show_plot:
        # display the initial image
        image.show(20)

    print('computing initial objective function value...')
    print('objective function value: %e' % (obj_fun.value(image)))

    if verbose:
        disp = 3
        if opt:
            print('NOTE: below f(x) is the negative of the objective function value')
    else:
        disp = 0
    eps = 1e-6 # single precision round-off error level

    for iter in range(steps):

        # obtain gradient for subset = iter
        grad = obj_fun.get_subset_gradient(image, iter % 12)
        # zero the gradient outside the cylindric FOV
        filter.apply(grad)

        # compute step size bazed on an estimate of the largest
        # eigenvalue lmd_max of the Hessian H
        # note that lmd_max = max |H v|/|v|
        if iter == 0:
            image0 = image
            grad0 = grad
            # in the quadratic case F(v) = (H v, v)/2,
            # grad F(v) = H v, hence a rough idea about lmd_max 
            # is given by
            lmd_max = 2*grad.norm()/image.norm()
            tau = 1/lmd_max
            maxstep = tau
        else:
            di = image - image0
            dg = grad - grad0 
            # dg = H di, hence a rough idea about lmd_max is given by
            lmd_max = 2*dg.norm()/di.norm()
            # alternative smaller estimate for lmd_max is
            #lmd_max = -2*dg.dot(di)/di.dot(di)
            tau = min(maxstep, 1/lmd_max)

        if opt:
            # find the optimal step size tau
            fun = lambda x: -obj_fun.value(image + x*grad)
            tau = scipy.optimize.fminbound \
                (fun, 0, 2*maxstep, xtol = 1e-4, maxfun = 4, disp = disp)

        print('using step size %f' % tau)

        # perform truncated steepest descent step
        new_image = trunc(image + tau*grad)
        diff = new_image - image
        rc = diff.norm()/image.norm()
        print('step %d, change in image %e' % (iter, rc))
        image = new_image
        # filter the new image
        filter.apply(image)

        if show_plot:
            # display the current image estimate
            image.show(20)

    if not opt or disp == 0:
        print('computing attained objective function value...')
        print('objective function value: %e' % (obj_fun.value(image)))


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
