'''OSSPS reconstruction demo with user-controlled iterations

Usage:
  ossps_reconstruction [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -g <file>, --init=<file>    initial image guess file
                              [default: test_image_PM_QP_6.hv]
  -f <fact>, --penf=<fact>    penalty factor [default: 0]
  -s <subs>, --subs=<subs>    number of subsets [default: 4]
  -i <iter>, --iter=<iter>    number of iterations [default: 2]
  -e <engn>, --engn=<engn>    reconstruction engine [default: Stir]
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
exec('from p' + args['--engn'] + ' import *')

# process command-line options
pen_factor = args['--penf']
num_subsets = int(args['--subs'])
num_subiterations = int(args['--iter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)
init_file = args['--init']

def main():

    # output goes to files
    printer = Printer('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    am = AcquisitionModelUsingMatrix()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    ad = AcquisitionData(raw_data_file)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = make_Poisson_loglikelihood(ad)
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)
    fact = float(pen_factor)
    obj_fun.set_prior(QuadraticPrior().set_penalisation_factor(fact))

    # select Ordered Subsets Separable Paraboloidal Surrogate
    # as the reconstruction algorithm
    recon = OSSPSReconstructor()
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)
    recon.set_objective_function(obj_fun)

    # read an initial estimate for the reconstructed image from the file
    # specified by --init option
    init_image_file = existing_filepath(data_path, init_file)
    image = ImageData(init_image_file)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iter in range(num_subiterations):
        print('\n------------- Subiteration %d' % recon.get_subiteration_num())
        # perform an iteration
        recon.update(image)
        # plot the current image
        data = image.as_array()
        pylab.figure(iter + 1)
        pylab.imshow(data[10,:,:])
        print('close Figure %d window to continue' % (iter + 1))
        pylab.show()

    # display image
    image.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
