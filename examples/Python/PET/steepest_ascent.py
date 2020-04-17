'''Steepest ascent demo.
Applies few steps of steepest ascent for the maximization of Poisson log-likelihood
objective function using subset gradients.

Usage:
  steepest_ascent [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -x <step>, --step=<step>     steepest ascent step size parameter,
                              use a negative value to opt for the optimal value
                              [default: 0.5]
  -s <nstp>, --steps=<nstp>   number of steepest descent steps [default: 3]
  -v, --verbose               verbose
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
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

import scipy.optimize

from pUtilities import show_2D_array

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

# process command-line options
step = float(args['--step'])
steps = int(args['--steps'])
verbose = args['--verbose']
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)

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

    # display the initial image
    image_as_3D_array = image.as_array()
    show_2D_array('Initial image', image_as_3D_array[20,:,:])

    print('computing initial objective function value...')
    print('objective function value: %e' % (obj_fun.value(image)))

    if verbose:
        disp = 3
        if step < 0:
            print('NOTE: below f(x) is the negative of the objective function value')
    else:
        disp = 0
    eps = 1e-6 # single precision round-off error level

    for iter in range(steps):

        # obtain gradient for subset = iter
        grad = obj_fun.get_subset_gradient(image, iter)
        # zero the gradient outside the cylindric FOV
        filter.apply(grad)
        grad_as_3D_array = grad.as_array()

        max_image = image_as_3D_array.max()
        max_grad = abs(grad_as_3D_array).max()
        delta = max_grad*eps

        # find maximal steepest ascent step parameter x in image + x*grad 
        # such that the new image remains positive;
        # since image is non-negative, the step size is limited by negative
        # gradients: it should not exceed -image/grad = abs(image/grad) at
        # points where grad is negative, thus, maximal x is min(abs(image/grad))
        # taken over such points

        # avoid division by zero at the next step
        grad_as_3D_array[abs(grad_as_3D_array) < delta] = delta
        # take the ratio of image to gradient
        ratio = abs(image_as_3D_array/grad_as_3D_array)
        # select points that are (i) inside cylindric FOV and (ii) such that
        # the gradient at them is negative
        select = numpy.logical_and(image_as_3D_array > 0, grad_as_3D_array < 0)
        if select.any():
            # at least one point is selected
            maxstep = ratio[select].min()
        else:
            # no such points - use a plausible value based on 'step' and
            # image-to-gradient ratio
            maxstep = abs(step)*max_image/max_grad

        # at some voxels image values may be close to zero and the gradient may
        # also be close to zero there; hence, both may become negative because
        # of round-off errors;
        # find such voxels and freeze them
        exclude = numpy.logical_and(image_as_3D_array <= 0, grad_as_3D_array < 0)
        grad_as_3D_array[exclude] = 0
        grad.fill(grad_as_3D_array)

        if step < 0:
            # find the optimal step size x
            fun = lambda x: -obj_fun.value(image + x*grad)
            x = scipy.optimize.fminbound \
                (fun, 0, maxstep, xtol = 1e-4, maxfun = 3, disp = disp)
        else:
            # x is such that the relative change in image is not greater than 'step'
            x = step*max_image/max_grad
            if x > maxstep:
                x = maxstep

        # perform steepest descent step
        print('step %d, max change in image %e' % (iter, x*max_grad))
        image = image + x*grad
        # filter the new image
        filter.apply(image)

        # display the current image estimate
        image_as_3D_array = image.as_array()
        show_2D_array('Current image', image_as_3D_array[20,:,:])

        # quit if the new image has substantially negative values
        min_image = image_as_3D_array.min()
        if min_image < -eps:
            print('image minimum is negative: %e' % min_image)
            break

    if step > 0 or disp == 0:
        print('computing attained objective function value...')
        print('objective function value: %e' % (obj_fun.value(image)))

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
