'''User-implemented OSMAPOSL reconstruction

Usage:
  user_osmaposl [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -s <subs>, --subs=<subs>    number of subsets [default: 12]
  -i <siter>, --subiter=<siter>    number of sub-iterations [default: 2]
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
exec('from p' + args['--engine'] + ' import *')

# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)

# user implementation of Ordered Subset Maximum A Posteriori One Step Late
# reconstruction algorithm
def my_osmaposl(image, obj_fun, prior, filter, num_subsets, num_subiterations):

    for sub_iter in range(1, num_subiterations + 1):
        print('\n------------- Subiteration %d' % sub_iter) 

        # select subset
        subset = (sub_iter - 1) % num_subsets

        # get sensitivity as ImageData
        sens_image = obj_fun.get_subset_sensitivity(subset)

        # get backprojection of the ratio of measured to estimated
        # acquisition data)
        grad_image = obj_fun.get_backprojection_of_acquisition_ratio\
                     (image, subset)

        # get gradient of prior as ImageData
        prior_grad_image = prior.get_gradient(image)

        # update image data
        denom = sens_image + prior_grad_image/num_subsets
        update = grad_image/denom
        image = image*update

        # apply filter
        filter.apply(image)

    return image

def main():

    # output goes to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # create acquisition model
    acq_model = AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # create filter that zeroes the image outside a cylinder of the same
    # diameter as the image xy-section size
    filter = TruncateToCylinderProcessor()

    # create initial image estimate
    image_size = (111, 111, 31)
    voxel_size = (3, 3, 3.375) # voxel sizes are in mm
    image = ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)

    # create prior
    prior = QuadraticPrior()
    prior.set_penalisation_factor(0.5)
    prior.set_up(image)

    # create objective function
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_num_subsets(num_subsets)
    obj_fun.set_up(image)

    # reconstruct using your own SIRF-based implementation of OSMAPOSL
    image = my_osmaposl \
        (image, obj_fun, prior, filter, num_subsets, num_subiterations)

    # show reconstructed image at z = 20
    image_array = image.as_array()
    show_2D_array('Reconstructed image at z = 20', image_array[20,:,:])

    image.write('my_image.hv')

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
