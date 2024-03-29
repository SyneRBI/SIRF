'''User-implemented OSMAPOSL reconstruction

Usage:
  user_osmaposl [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -s <subs>, --subs=<subs>    number of subsets [default: 12]
  -i <sit>, --subiter=<sit>   number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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

from sirf.Utilities import error, examples_data_path, existing_filepath
from sirf.Utilities import show_2D_array

# import engine module
import importlib
engine = args['--engine']
pet = importlib.import_module('sirf.' + engine)


# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
show_plot = not args['--non-interactive']


# user implementation of Ordered Subset Maximum A Posteriori One Step Late
# reconstruction algorithm
def my_osmaposl(image, obj_fun, prior, Filter, num_subsets, num_subiterations):

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
        Filter.apply(image)

    return image


def main():

    # output goes to files
    _ = pet.MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # create acquisition model
    acq_model = pet.AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = pet.AcquisitionData(raw_data_file)

    # create filter that zeroes the image outside a cylinder of the same
    # diameter as the image xy-section size
    Filter = pet.TruncateToCylinderProcessor()

    # create initial image estimate
    image_size = (31, 111, 111)
    voxel_size = (3.375, 3, 3) # voxel sizes are in mm
    image = pet.ImageData()
    image.initialise(image_size, voxel_size)
    image.fill(1.0)

    # create prior
    prior = pet.QuadraticPrior()
    prior.set_penalisation_factor(0.5)
    prior.set_up(image)

    # create objective function
    obj_fun = pet.make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_num_subsets(num_subsets)
    obj_fun.set_up(image)

    # reconstruct using your own SIRF-based implementation of OSMAPOSL
    image = my_osmaposl \
        (image, obj_fun, prior, Filter, num_subsets, num_subiterations)

    if show_plot:
        # show reconstructed image at z = 20
        image_array = image.as_array()
        show_2D_array('Reconstructed image at z = 20', image_array[20,:,:])

#    image.write('my_image.hv')


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
