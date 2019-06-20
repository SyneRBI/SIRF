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

    # Run the scatter simulation
    sss.run_scatter_simulation()

    # Get the simulated sinogram.
    # Note: The sinogram has been downscaled.
    sss_data = sss.get_simulated_data()

    # show simulated scatter data
    simulated_scatter_as_array = sss_data.as_array()
    show_2D_array('Scatter simulation', simulated_scatter_as_array[0,0,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)