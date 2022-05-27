'''

Forward projection demo: creates an image, projects it to simulate
acquisition data and backprojects

Usage:
  acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: template_sinogram.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/SPECT
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

## CCP SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017, 2019, 2022 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (SyneRBI)
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
import sirf.STIR

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = sirf.STIR.examples_data_path('SPECT')
raw_data_file = sirf.STIR.existing_filepath(data_path, data_file)
output_file = args['--output']

def create_sample_image(image):
    '''fill the image with some simple geometric shapes.'''
    image.fill(0)
    # create a shape
    shape = sirf.STIR.EllipticCylinder()
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

def main():

    ## AcquisitionData.set_storage_scheme('mem')

    # no info printing from the engine, warnings and errors sent to stdout
    # msg_red = MessageRedirector()
    # output goes to files
    msg_red = sirf.STIR.MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # raw data to be used as a template for the acquisition model
    acq_template = sirf.STIR.AcquisitionData(raw_data_file)

    # create image with suitable sizes
    image = acq_template.create_uniform_image()
    create_sample_image(image)
    image.write("simulated_image.hv")

    # show the phantom image
    image_array = image.as_array()
    show_2D_array('Phantom image', image_array[0,:,:])

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model_matrix = sirf.STIR.SPECTUBMatrix();
    acq_model = sirf.STIR.AcquisitionModelUsingMatrix(acq_model_matrix)

    # require same number slices and equal z-sampling for projection data & image
    image = image.zoom_image(zooms=(0.5, 1.0, 1.0))
    print('projecting image...')
    # project the image to obtain simulated acquisition data
    # data from raw_data_file is used as a template
    acq_model.set_up(acq_template, image)
    simulated_data = acq_template.get_uniform_copy()
    acq_model.forward(image, 0, 1, simulated_data)
    # simulated_data = acq_model.forward(image, 0, 4)
    if output_file is not None:
        simulated_data.write(output_file)

    # show simulated acquisition data
    simulated_data_as_array = simulated_data.as_array()
    show_2D_array('Forward projection', simulated_data_as_array[0, 0,:,:])

    print('backprojecting the forward projection...')
    # backproject the computed forward projection
    back_projected_image = acq_model.backward(simulated_data, 0, 1)

    back_projected_image_as_array = back_projected_image.as_array()
    show_2D_array('Backprojection', back_projected_image_as_array[0,:,:])

try:
    main()
    print('done')
except sirf.STIR.error as err:
    print('%s' % err.value)
