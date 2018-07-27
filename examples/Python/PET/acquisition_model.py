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

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2018 University College London.
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
addv = float(args['--addv'])
back = float(args['--back'])
beff = 1/float(args['--norm'])
output_file = args['--output']

def main():

##    AcquisitionData.set_storage_scheme('mem')

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = MessageRedirector()
    # output goes to files
##    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # create an empty image
    image = ImageData()
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375) # voxel sizes are in mm
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

    # z-pixel coordinate of the xy-crossection to show
    z = int(image_size[2]/2)

    # show the phantom image
    image_array = image.as_array()
    show_2D_array('Phantom image', image_array[z,:,:])

    # raw data to be used as a template for the acquisition model
    acq_template = AcquisitionData(raw_data_file)

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix()

    # testing bin efficiencies
    bin_eff = acq_template.clone()
    bin_eff.fill(beff)
    bin_eff_arr = bin_eff.as_array()
    # As an example, if bin efficiencies are non-trivial, set a portion of them to zero;
    # this should zero the corresponding portion of forward projection
    # and 'damage' the backprojection making it look less like the
    # actual image
    if beff != 1:
        bin_eff_arr[:,10:50,:] = 0
    show_2D_array('Bin efficiencies', bin_eff_arr[0,:,:])
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
    acq_model.forward(image, 0, 4, simulated_data)
#    simulated_data = acq_model.forward(image, 0, 4)
    if output_file is not None:
        simulated_data.write(output_file)

    # show simulated acquisition data
    simulated_data_as_array = simulated_data.as_array()
    show_2D_array('Forward projection', simulated_data_as_array[0,:,:])

    print('backprojecting the forward projection...')
    # backproject the computed forward projection
    # note that the backprojection takes the acquisition sensitivy model asm into account as well
    back_projected_image = acq_model.backward(simulated_data, 0, 4)

    back_projected_image_as_array = back_projected_image.as_array()
    show_2D_array('Backprojection', back_projected_image_as_array[z,:,:])
    
    acq_model.direct(image, 0, 4, simulated_data)
#    simulated_data = acq_model.forward(image, 0, 4)
    if output_file is not None:
        simulated_data.write(output_file)

    # show simulated acquisition data
    simulated_data_as_array_direct = simulated_data.as_array()
    show_2D_array('Direct projection', simulated_data_as_array_direct[0,:,:])
    back_projected_image_adj = acq_model.adjoint(simulated_data, 0, 4)

    back_projected_image_as_array_adj = back_projected_image_adj.as_array()
    show_2D_array('Adjoint projection', back_projected_image_as_array_adj[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
