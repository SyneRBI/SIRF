'''OSSPS reconstruction demo

Usage:
  ossps_reconstruction [--help | options]

Options:
  -d <file>, --file=<file>     raw data file
                               [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
  -g <file>, --init=<file>     initial image guess file
                               [default: test_image_PM_QP_6.hv]
  -f <fact>, --penf=<fact>     penalty factor [default: 1000]
  -s <subs>, --subs=<subs>     number of subsets [default: 4]
  -i <sitr>, --subiter=<sitr>  number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  --non-interactive            do not show plots
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

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')


# process command-line options
pen_factor = args['--penf']
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
init_file = args['--init']
show_plot = not args['--non-interactive']


def main():

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = MessageRedirector()
    # all engine's printing goes to files
##    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # read an initial estimate for the reconstructed image from the file
    # specified by --init option
    init_image_file = existing_filepath(data_path, init_file)
    image = ImageData(init_image_file)
    image_array = image.as_array()
    if show_plot:
        show_2D_array('Initial image', image_array[10,:,:])

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    fact = float(pen_factor)
    obj_fun.set_prior(QuadraticPrior().set_penalisation_factor(fact))

    # select Ordered Subsets Separable Paraboloidal Surrogate
    # as the reconstruction algorithm
    recon = OSSPSReconstructor()
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)
    recon.set_objective_function(obj_fun)
    recon.set_input(acq_data)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # open up the user's access to  the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iter in range(num_subiterations):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform a sub-iteration
        recon.update(image)
        if not show_plot:
            continue
        # display the current image at z = 10
        image_array = image.as_array()
        pylab.figure(iter + 1)
        pylab.imshow(image_array[10,:,:])
        print('You may neet to close Figure %d window to continue' % (iter + 1))
    if show_plot:
        pylab.show()

    if show_plot:
        # display the reconstructed image
        #image.show(15) # show one
        #image.show([1, 10, 11, 12]) # show list
        #image.show((1, 10, 11, 12))
        #image.show(range(2, 30, 3)) # show range
        image.show(title='Reconstructed images') # show all


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
