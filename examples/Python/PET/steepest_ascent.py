'''Steepest ascent demo

Usage:
  steepest_ascent [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -t <step>, --tau=<step>     steepest ascent step size parameter,
                              use a negative value to opt for the optimal value
                              [default: -1]
  -s <nstp>, --steps=<nstp>   number of steepest descent steps [default: 3]
  -v, --verbose               verbose
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
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

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
tau = float(args['--tau'])
steps = int(args['--steps'])
verbose = args['--verbose']
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)

def show(fig, title, data):
    pylab.figure(fig)
    pylab.title(title)
    pylab.imshow(data)
    pylab.colorbar()
    print('close window to continue')
    pylab.show()

def main():

    # output goes to files
    printer = Printer('info.txt', 'warn.txt', 'errr.txt')

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
    image_size = (nx, ny, nz)
    voxel_size = (3, 3, 3.375) # sizes are in mm
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)
    # apply the filter to the image
    filter.apply(image)

    # create objective function of Poisson logarithmic likelihood type
    # compatible with the acquisition data type
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_acquisition_data(acq_data)
    obj_fun.set_num_subsets(12)
    obj_fun.set_up(image)

    # display the initial image
    image_as_3D_array = image.as_array()
    show(1, 'initial image', image_as_3D_array[20,:,:])

    print('computing initial objective function value...')
    print('objective function value: %e' % (obj_fun.value(image)))

    if verbose:
        disp = 3
        print('NOTE: below f(x) is the negative of the objective function value')
    else:
        disp = 0
    eps = 1e-6

    for iter in range(steps):

        # obtain gradient
        grad = obj_fun.gradient(image, 0)
        filter.apply(grad)
        grad_as_3D_array = grad.as_array()

        # find maximal steepest descent step parameter t in image + t*grad 
        # such that the new image remains positive
        max_image = image_as_3D_array.max()
        max_grad = abs(grad_as_3D_array).max()
        grad_as_3D_array[abs(grad_as_3D_array) < eps] = eps
        r = image_as_3D_array/grad_as_3D_array
        grad_as_3D_array[image_as_3D_array <= 0] = 0
        d = r[r < 0]
        if d.shape[0] > 0:
            maxstep = -d.max()
        else:
            maxstep = abs(tau)*max_image/max_grad

        if tau < 0:
            # find the optimal t
            fun = lambda t: -obj_fun.value \
                  (image.fill(image_as_3D_array + t*grad_as_3D_array))
            t = scipy.optimize.fminbound \
                (fun, 0, maxstep, xtol = 1e-4, maxfun = 3, disp = disp)
        else:
            # t is such that the relative change in image is not greater than tau
            t = tau*max_image/max_grad
            if t > maxstep:
                t = maxstep

        # perform steepest descent step
        print('step %d, max change in image %e' % (iter, t*max_grad))
        image_as_3D_array = image_as_3D_array + t*grad_as_3D_array

        # filter the new image
        image.fill(image_as_3D_array)
        filter.apply(image)
        image_as_3D_array = image.as_array()

        # display the current image estimate
        show(iter + 1, 'current image', image_as_3D_array[20,:,:])

        # quit if the new image has negative values
        min_image = image_as_3D_array.min()
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
