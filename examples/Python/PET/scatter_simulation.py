'''Scatter simulation demo: creates an image, projects it to simulate
acquisition data and runs a scatter simulation

Usage:
  scatter_simulation [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: scatter_template.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data [default: scatter_output.hs]

WARNING: Currently this computes scatter at the same sampling as the template. If you use
real projection data as template, this will be very slow (and might run out of memory).
(The STIR upsampling facility is currently not yet in SIRF. Sorry)
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2019 University of Hull
## Copyright 2020 University College London
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
import matplotlib.pyplot as plt

args = docopt(__doc__, version=__version__)

# import engine module
import sirf.STIR as PET
import os
from sirf.Utilities import show_2D_array

# process command-line options
data_file = args['--file']
data_path = args['--path']
print ("FILE",__file__)
if data_path is None:
    data_path = os.path.join(os.path.dirname(__file__), '..', '..', 'parameter_files')
    print(data_path)
acq_template_filename = PET.existing_filepath(data_path, data_file)
output_file = args['--output']


def main():
    ##    PET.AcquisitionData.set_storage_scheme('memory')

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = PET.MessageRedirector()

    # Create a template Acquisition Model
    #acq_template = AcquisitionData('Siemens mMR', 1, 0, 1)
    acq_template = PET.AcquisitionData(acq_template_filename)#q.get_uniform_copy()

    # create the attenuation image
    atten_image = PET.ImageData(acq_template)
    image_size = atten_image.dimensions()
    voxel_size = atten_image.voxel_sizes()

    # create a cylindrical water phantom
    water_cyl = PET.EllipticCylinder()
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
    act_image = atten_image.clone()
    act_image.fill(0.0)

    # create the activity cylinder
    act_cyl = PET.EllipticCylinder()
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
    sss = PET.SingleScatterSimulator()

    # Set the attenuation image
    sss.set_attenuation_image(atten_image)

    # set-up the scatter simulator
    sss.set_up(acq_template, act_image)
    # Simulate!
    sss_data = sss.forward(act_image)

    # show simulated scatter data
    simulated_scatter_as_array = sss_data.as_array()
    show_2D_array('scatter simulation', simulated_scatter_as_array[0,0,:,:])

    sss_data.write(output_file)

    ## let's also compute the unscattered counts (at the same low resolution) and compare

    acq_model = PET.AcquisitionModelUsingRayTracingMatrix()
    asm = PET.AcquisitionSensitivityModel(atten_image, acq_model)
    acq_model.set_acquisition_sensitivity(asm)

    acq_model.set_up(acq_template, act_image)
    #unscattered_data = acq_template.get_uniform_copy()
    unscattered_data = acq_model.forward(act_image)
    simulated_unscatter_as_array = unscattered_data.as_array()
    show_2D_array('unscattered simulation', simulated_unscatter_as_array[0,0,:,:])

    plt.figure()
    ax = plt.subplot(111)
    plt.plot(simulated_unscatter_as_array[0,4,0,:], label='unscattered')
    plt.plot(simulated_scatter_as_array[0,4,0,:], label='scattered')
    ax.legend()
    plt.show()
try:
    main()
    print('done')
except PET.error as err:
    print('%s' % err.value)
