'''Forward projection demo: creates an image, projects it to simulate
acquisition data and backprojects

Usage:
  acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2019 University College London.
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

from pUtilities import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)
output_file = args['--output']

def main():

##    AcquisitionData.set_storage_scheme('mem')

    # no info printing from the engine, warnings and errors sent to stdout
#    msg_red = MessageRedirector()
    # output goes to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # create an empty image
    image = ImageData()
    image_size = (111, 111, 64)
    voxel_size = (3, 3, 3.32) # voxel sizes are in mm
    image.initialise(image_size, voxel_size)

    # create a shape
    shape = EllipticCylinder()
    shape.set_length(400)
    shape.set_radii((100, 40))
    shape.set_origin((0, 60, 10))

    # add the shape to the image
    image.add_shape(shape, scale = 1)

    # add another shape
    shape.set_radii((30, 30))
    shape.set_origin((60, -30, 10))
    image.add_shape(shape, scale = 1.5)

    # add another shape
    shape.set_origin((-60, -30, 10))
    image.add_shape(shape, scale = 0.75)

    image.write("simulated_image.hv")
    # z-pixel coordinate of the xy-cross-section to show
    z = int(image_size[2]/2)

    # show the phantom image
    image_array = image.as_array()
    show_2D_array('Phantom image', image_array[z,:,:])

    # raw data to be used as a template for the acquisition model
    acq_template = AcquisitionData(raw_data_file)

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model_matrix = SPECTUBMatrix();
    acq_model = AcquisitionModelUsingMatrix(acq_model_matrix)

    print('projecting image...')
    # project the image to obtain simulated acquisition data
    # data from raw_data_file is used as a template
    acq_model.set_up(acq_template, image)
    simulated_data = acq_template.get_uniform_copy()
    acq_model.forward(image, 0, 1, simulated_data)
#    simulated_data = acq_model.forward(image, 0, 4)
    if output_file is not None:
        simulated_data.write(output_file)

    # show simulated acquisition data
    simulated_data_as_array = simulated_data.as_array()
    middle_slice=simulated_data_as_array.shape[0]//2
    show_2D_array('Forward projection', simulated_data_as_array[middle_slice,:,:])

    print('backprojecting the forward projection...')
    # backproject the computed forward projection
    back_projected_image = acq_model.backward(simulated_data, 0, 1)

    back_projected_image_as_array = back_projected_image.as_array()
    show_2D_array('Backprojection', back_projected_image_as_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
