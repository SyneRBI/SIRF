'''Scatter simulation demo: creates an image, projects it to simulate
acquisition data and runs a scatter simulation

Usage:
  scatter_simulation [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -a <addv>, --addv=<addv>    additive term value [default: 0]
  -b <back>, --back=<back>    background term value [default: 0]
  -n <norm>, --norm=<norm>    normalization value [default: 1]
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2019 University of Hull
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

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

from pUtilities import show_2D_array

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
beff = 1 / float(args['--norm'])
output_file = args['--output']


def main():
    ##    AcquisitionData.set_storage_scheme('mem')

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = MessageRedirector()
    # output goes to files
    ##    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # Create a template Acquisition Model
    tmpl_acq_data = AcquisitionData('ECAT 931', 1, 0, 1)

    # create the attenuation image
    atten_image = ImageData(tmpl_acq_data)
    image_size = atten_image.dimensions()
    voxel_size = atten_image.voxel_sizes()

    # create a cylindrical water phantom
    water_cyl = EllipticCylinder()
    water_cyl.set_length(image_size[0]*voxel_size[0])
    water_cyl.set_radii((image_size[1]*voxel_size[1]*0.25, \
                     image_size[2]*voxel_size[2]*0.25))
    water_cyl.set_origin((image_size[0]*voxel_size[0]*0.5, 0, 0))

    # add the shape to the image
    atten_image.add_shape(water_cyl, scale = 9.687E-02)

    # z-pixel coordinate of the xy-crossection to show
    z = int(image_size[0]*0.5)

    # show the phantom image
    atten_image_array = atten_image.as_array()
    show_2D_array('Attenuation image', atten_image_array[z,:,:])

    # Create the activity image
    act_image = atten_image.clone();
    act_image.fill(0.0)

    # create the activity cylinder
    act_cyl = EllipticCylinder()
    act_cyl.set_length(image_size[0] * voxel_size[0])
    act_cyl.set_radii((image_size[1] * voxel_size[1] * 0.125, \
                         image_size[2] * voxel_size[2] * 0.125))
    act_cyl.set_origin((0, image_size[1] * voxel_size[1] * 0.06, \
                          image_size[2] * voxel_size[2] * 0.06))

    # add the shape to the image
    act_image.add_shape(act_cyl, scale=1)

    # z-pixel coordinate of the xy-crossection to show
    z = int(image_size[0] * 0.5)

    # show the phantom image
    act_image_array = act_image.as_array()
    show_2D_array('Activity image', act_image_array[z, :, :])

    # Create the Single Scatter Simulation model
    sss = SingleScatterSimulator()

    # Input in the Simulator information about the Acquisition
    sss.set_acquisition_data(tmpl_acq_data)

    # Set the attenuation image
    sss.set_attenuation_image(atten_image)

    # Set the activity image
    sss.set_activity_image(act_image)

    sss.run_scatter_simulation()

    # sss_data = sss.get_output()

    debug_stop = 0


#     # select acquisition model that implements the geometric
#     # forward projection by a ray tracing matrix multiplication
#     acq_model = AcquisitionModelUsingRayTracingMatrix()
#
#     # testing bin efficiencies
#     bin_eff = acq_template.clone()
#     bin_eff.fill(beff)
#     bin_eff_arr = bin_eff.as_array()
#     # As an example, if bin efficiencies are non-trivial, set a portion of them to zero;
#     # this should zero the corresponding portion of forward projection
#     # and 'damage' the backprojection making it look less like the
#     # actual image
#     if beff != 1:
#         bin_eff_arr[0,:,10:50,:] = 0
#     show_2D_array('Bin efficiencies', bin_eff_arr[0,0,:,:])
#     bin_eff.fill(bin_eff_arr)
#
#     asm = AcquisitionSensitivityModel(bin_eff)
#     acq_model.set_acquisition_sensitivity(asm)
#
#     # As an example, add both an additive term and background term
#     # (you normally wouldn't do this for real data)
#     add = acq_template.clone()
#     add.fill(addv)
#     acq_model.set_additive_term(add)
#
#     bck = acq_template.clone()
#     bck.fill(back)
#     acq_model.set_background_term(bck)
#
#     print('projecting image...')
#     # project the image to obtain simulated acquisition data
#     # data from raw_data_file is used as a template
#     acq_model.set_up(acq_template, image)
#     simulated_data = acq_template.get_uniform_copy()
#     acq_model.forward(image, 0, 4, simulated_data)
# #    simulated_data = acq_model.forward(image, 0, 4)
#     if output_file is not None:
#         simulated_data.write(output_file)
#
#     # show simulated acquisition data
#     simulated_data_as_array = simulated_data.as_array()
#     show_2D_array('Forward projection', simulated_data_as_array[0,0,:,:])
#
#     print('backprojecting the forward projection...')
#     # backproject the computed forward projection
#     # note that the backprojection takes the acquisition sensitivy model asm into account as well
#     back_projected_image = acq_model.backward(simulated_data, 0, 4)
#
#     back_projected_image_as_array = back_projected_image.as_array()
#     show_2D_array('Backprojection', back_projected_image_as_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)