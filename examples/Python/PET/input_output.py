'''PET I/O demo

Usage:
  input_output [--help | options]

Options:
  -p <path> , --path=<path>   path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -t <file> , --tfile=<file>  termplate file [default: my_forward_projection.hs]
  -a <name> , --afile=<name>  file to store simulated acquisition data
                              [default: simulated_data]
  -i <name> , --ifile=<name>  file to store phantom image data [default: phantom]
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

from pUtilities import show_2D_array

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
templ_file = args['--tfile']
templ_file = existing_filepath(data_path, templ_file)
acq_file = args['--afile']
img_file = args['--ifile']

def main():

    # engine's messages go to files, except error messages, which go to stdout
    msg_red = MessageRedirector('info.txt', 'warn.txt')

    acq_template = AcquisitionData(templ_file)
##    # create acquisition data from scanner parameters to be used as a template
##    print('creating Siemens_mMR acquisition data...')
##    acq_template = AcquisitionData('Siemens_mMR')
##    acq_dim = acq_template.dimensions()
##    print('acquisition data dimensions: %d sinograms, %d views, %d tang. pos.' \
##          % acq_dim)
##    # rebin to reduce the acquisition data size
##    print('rebinning...')
##    acq_template = acq_template.rebin(15)
    acq_dim = acq_template.dimensions()
    print('acquisition data dimensions: ' + \
          '%d TOF bins %d (non-TOF) sinograms, %d views, %d tang. pos.' \
          % acq_dim)

    # create image of dimensions and voxel sizes compatible with the scanner
    # geometry (stored in the AcquisitionData object ad)
    # and initialize each voxel to 1.0
    print('creating compatible phantom')
    image = acq_template.create_uniform_image()
    # show the image
    nz, ny, nx = image.dimensions()
    vz, vy, vx = image.voxel_sizes()
    print('phantom dimensions: %dx%dx%d' % (nz, ny, nx))
    print('phantom voxel sizes: %fx%fx%f' % (vz, vy, vx))
    image_size = (int(nz), 111, 111)
    voxel_size = (float(vz), 3.0, 3.0)
    image = ImageData()
    image.initialise(image_size, voxel_size)

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

    image_array = image.as_array()
    z = int(image_array.shape[0]/2)
    show_2D_array('Phantom', image_array[z,:,:])

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix()
    print('setting up the acquisition model...')
    acq_model.set_up(acq_template, image)
    # project the image to obtain simulated acquisition data
    print('projecting the phantom to create simulated acquisition data...')
    simulated_data = acq_model.forward(image)

    # copy the acquisition data into a Python array
    acq_array = simulated_data.as_array()
    acq_dim = acq_array.shape
##    print('acquisition data dimensions: %dx%dx%d' % acq_dim)
    z = acq_dim[1]//2
    show_2D_array('Simulated acquisition data', acq_array[0,z,:,:])

    # write acquisition data and image to files
    print('writing acquisition data...')
    simulated_data.write(acq_file)
    print('writing image...')
    image.write(img_file)

    # read acquisition data and image from files
    acq_data = AcquisitionData('simulated_data.hs')
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
##    print('acquisition data dimensions: %dx%dx%d' % acq_dim)
    z = acq_dim[1]//2
    show_2D_array('Simulated acquisition data', acq_array[0,z,:,:])

    # show the image again
    img = ImageData()
    img.read_from_file('phantom.hv')
    image_array = img.as_array()
##    print('phantom dimensions: %dx%dx%d' % image_array.shape[2::-1])
    z = int(image_array.shape[0]/2)
    show_2D_array('Phantom', image_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
