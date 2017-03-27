'''Forward projection demo: creates an image, projects it to simulate
acquisition data and displays

Usage:
  using_acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
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

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
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

    # create an empty image
    image = ImageData()
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375) # voxel sizes are in mm
    image.initialise(image_size, voxel_size)

    # create a shape
    shape = EllipticCylinder()

    # add a shape to the image
    shape.set_length(400)
    shape.set_radii((100, 40))
    shape.set_origin((0, 60, 10))
    image.add_shape(shape, scale = 1)

    # add another shape
    shape.set_radii((30, 30))
    shape.set_origin((60, -30, 10))
    image.add_shape(shape, scale = 1.5)

    # add another shape
    shape.set_origin((-60, -30, 10))
    image.add_shape(shape, scale = 0.75)

    # z-pixel coordinate of the xy-crossection to show
    z = int(image_size[2]/2)

    # show the phantom image
    data = image.as_array()
    show(1, 'Figure 1: phantom image', data[z,:,:])

    # select acquisition model that implements the geometric
    # forward projection by a matrix multiplication;
    # matrix type defaults to ray tracing
    acq_model = AcquisitionModelUsingMatrix()

    print('projecting image...')
    # project the image to obtain simulated acquisition data
    # data from raw_data_file is used as a template
    template = AcquisitionData(raw_data_file)
    acq_model.set_up(template, image)
    simulated_data = acq_model.forward(image)
    # if the projection data is very large, it can be stored in a file
    # simulated_data = am.forward(image, 'proj_data.hs')

    # show simulated acquisition data
    simulated_data_as_array = simulated_data.as_array()
    show(2, 'Figure 2: forward projection', simulated_data_as_array[z,:,:])

    print('backprojecting the forward projection...')
    # backproject the computed forward projection
    back_projected_image = acq_model.backward(simulated_data)

    back_projected_image_as_array = back_projected_image.as_array()
    show(3, 'Figure 3: back projection', back_projected_image_as_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('exception occured: %s' % err.value)
