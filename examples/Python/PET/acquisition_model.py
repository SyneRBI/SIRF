'''Forward projection demo: creates an image, projects it to simulate
acquisition data and backprojects

Usage:
  acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -a <addv>, --addv=<addv>    additive term value [default: 0]
  -b <back>, --back=<back>    background term value [default: 0]
  -n <norm>, --norm=<norm>    normalization value [default: 1]
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
  --parallelproj              use Parallelproj for projections
  --non-interactive           do not show plots

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
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

from sirf.Utilities import show_2D_array

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')


# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
addv = float(args['--addv'])
back = float(args['--back'])
beff = 1/float(args['--norm'])
output_file = args['--output']
parallelproj = args['--parallelproj']
show_plot = not args['--non-interactive']


try:
    import matplotlib.pyplot as plt
except:
    show_plot = False


def main():

    print(scanner_names())

##    AcquisitionData.set_storage_scheme('mem')

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = MessageRedirector()
    # output goes to files
##    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # raw data to be used as a template for the acquisition model
    acq_template = AcquisitionData(raw_data_file)

    # create an empty image
    image = acq_template.create_uniform_image(0.0, xy=111)
    image_size = image.dimensions()
    print('image size: %d by %d by %d' % image_size)

    # create a shape
    shape = EllipticCylinder()
    shape.set_length(400)
    shape.set_radii((40, 100))
    shape.set_origin((10, 60, 0))

    # add the shape to the image
    image.add_shape(shape, scale = 1)

    # add another shape
    shape.set_radii((30, 30))
    shape.set_origin((10, -30, 60))
    image.add_shape(shape, scale = 1.5)

    # add another shape
    shape.set_origin((10, -30, -60))
    image.add_shape(shape, scale = 0.75)

    # apply Gaussian filter
    filter = SeparableGaussianImageFilter()
    filter.set_fwhms((10, 20, 30))
    filter.set_max_kernel_sizes((10, 10, 2))
    filter.set_normalise()
    filter.set_up(image)
    filter.apply(image)

    # z-pixel coordinate of the xy-crossection to show
    z = int(image_size[0]/2)

    if show_plot:
        # show the phantom image
        image_array = image.as_array()
        show_2D_array('Phantom image', image_array[z,:,:])

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    if parallelproj:
        acq_model = AcquisitionModelUsingParallelproj()
        num_subsets = 1
    else:
        acq_model = AcquisitionModelUsingRayTracingMatrix()
        num_subsets = 4

    # testing bin efficiencies
    bin_eff = acq_template.clone()
    bin_eff.fill(beff)
    bin_eff_arr = bin_eff.as_array()
    # As an example, if bin efficiencies are non-trivial, set a portion of them to zero;
    # this should zero the corresponding portion of forward projection
    # and 'damage' the backprojection making it look less like the
    # actual image
    if beff != 1:
        bin_eff_arr[0,:,10:50,:] = 0
    if show_plot:
        show_2D_array('Bin efficiencies', bin_eff_arr[0,0,:,:])
    bin_eff.fill(bin_eff_arr)

    asm = AcquisitionSensitivityModel(bin_eff)
    acq_model.set_acquisition_sensitivity(asm)

    # As an example, add both an additive term and background term
    # (you normally wouldn't do this for real data)
    add = acq_template.clone()
    add.fill(addv)
    acq_model.set_additive_term(add)

    bck = acq_template.clone()
    bck.fill(back)
    acq_model.set_background_term(bck)

    print('projecting image...')
    # project the image to obtain simulated acquisition data
    # data from raw_data_file is used as a template
    acq_model.set_up(acq_template, image)
    simulated_data = acq_template.get_uniform_copy()
    acq_model.forward(image, 0, num_subsets, simulated_data)
#    simulated_data = acq_model.forward(image, 0, num_subsets)
    if output_file is not None:
        simulated_data.write(output_file)

    print('\n--- Computing the norm of the linear part A of acquisition model...')
    acqm_norm = acq_model.norm()
    image_norm = image.norm()
    acqd_norm = simulated_data.norm()
    print('\n--- The computed norm is |A| = %f, checking...' % acqm_norm)
    print('    image data x norm: |x| = %f' % image_norm)
    print('    forward projected data A x norm: |A x| = %f' % acqd_norm)
    acqd_bound = acqm_norm*image_norm
    msg = '    |A x| must be less than or equal to |A||x| = %f'
    if acqd_norm <= acqd_bound:
        msg += ' - ok\n'
    else:
        msg += ' - ???\n'
    print( msg % acqd_bound)

    if show_plot:
        # show simulated acquisition data
        simulated_data_as_array = simulated_data.as_array()
        show_2D_array('Forward projection (subset 0)', simulated_data_as_array[0,0,:,:])

    print('backprojecting the forward projection...')
    # backproject the computed forward projection
    # note that the backprojection takes the acquisition sensitivy model asm into account as well
    back_projected_image = acq_model.backward(simulated_data, 0, num_subsets)
    back_projected_image_as_array = back_projected_image.as_array()
    if show_plot:
        show_2D_array('Backprojection', back_projected_image_as_array[z,:,:])

    # backproject again, this time into pre-allocated image
    back_projected_image.fill(0.0)
    acq_model.backward(simulated_data, 0, num_subsets, out=back_projected_image)
    back_projected_image_as_array = back_projected_image.as_array()
    if show_plot:
        msg = 'Backprojection into pre-allocated image'
        show_2D_array(msg, back_projected_image_as_array[z,:,:])

    # do same with pre-smoothing (often used for resolution modelling)
    print('Using some PSF modelling for comparison')
    smoother = SeparableGaussianImageFilter()
    smoother.set_fwhms((6,11,12))
    acq_model.set_image_data_processor(smoother)
    acq_model.set_up(acq_template, image)
    simulated_data_PSF = acq_template.get_uniform_copy()
    acq_model.forward(image, 0, num_subsets, simulated_data_PSF)
    if show_plot:
        simulated_data_PSF_as_array = simulated_data_PSF.as_array()
        plt.figure()
        plt.plot(simulated_data_as_array[0,0,0,:], label="no PSF")
        plt.plot(simulated_data_PSF_as_array[0,0,0,:], label="PSF")
        plt.title('Diff Forward projection without/ with smoothing (first view)')
        plt.legend()
    # backprojection
    back_projected_image_PSF = acq_model.backward(simulated_data, 0, num_subsets)
    if show_plot:
        back_projected_image_PSF_as_array = back_projected_image_PSF.as_array()
        y = back_projected_image_as_array.shape[1]//2;
        plt.figure()
        plt.plot(back_projected_image_as_array[z,y,:], label="no PSF")
        plt.plot(back_projected_image_PSF_as_array[z,y,:], label="PSF")
        plt.title('Diff Back projection without/ with smoothing (central horizontal line)')
        plt.legend()

    # direct is alias for the forward method for a linear AcquisitionModel
    # raises error if the AcquisitionModel is not linear.
    try:
        acq_model.num_subsets = num_subsets
        acq_model.direct(image, simulated_data)
    except error as err:
        print('%s' % err.value)
        print('Extracting the linear acquisition model...')
        lin_acq_model = acq_model.get_linear_acquisition_model()
        lin_acq_model.direct(image, simulated_data)

    if show_plot:
        # show simulated acquisition data
        simulated_data_as_array_direct = simulated_data.as_array()
        show_2D_array('Direct projection', simulated_data_as_array_direct[0,0,:,:])
    
    # adjoint is an alias for the backward method for a linear AcquisitionModel
    # raises error if the AcquisitionModel is not linear.
    try:
        back_projected_image_adj = acq_model.adjoint(simulated_data)
    except error as err:
        print('%s' % err.value)
        print('Extracting the linear acquisition model...')
        lin_acq_model = acq_model.get_linear_acquisition_model()
        back_projected_image_adj = lin_acq_model.adjoint(simulated_data)

    if show_plot:
        back_projected_image_as_array_adj = back_projected_image_adj.as_array()
        show_2D_array('Adjoint projection', back_projected_image_as_array_adj[z,:,:])

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    print('%s' % err.value)
