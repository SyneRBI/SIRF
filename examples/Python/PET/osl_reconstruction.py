'''One-Step-Late reconstruction demo.
We use the OSMAPOSL reconstructor in this demo. This reconstructor
implements an Ordered Subsets (OS) version of the One Step Late algorithm (OSL)
from Green et al for Maximum a Posteriori (MAP) maximisation. 
This is an algorithm often used in PET. However, it has a known 
problem: it can diverge if the penalty factor is too large.
(Try to run this with a small penalty and a very large one).

We illustrate OSL here with 2 priors: the quadratic prior and the
Parallel Level Sets prior. The latter uses an anatomical image to preserve
edges.

Usage:
  osl_reconstruction [--help | options]

Options:
  -d <file>, --file=<file>    raw data file [default: my_forward_projection.hs]
  -a <file>, --anim=<file>    anatomical image file
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -f <fact>, --penf=<fact>    penalty factor [default: 10]
  -s <subs>, --subs=<subs>    number of subsets [default: 12]
  -i <iter>, --subiter=<iter>    number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2018 University College London.
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

# import engine module
import importlib
engine = args['--engine']
pet = importlib.import_module('sirf.' + engine)


# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
pen_factor = args['--penf']
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
if args['--anim'] is not None:
    ai_file = existing_filepath(data_path, args['--anim'])
else:
    ai_file = None
show_plot = not args['--non-interactive']


def main():
 
    # direct all engine's information and warnings printing to files
    _ = pet.MessageRedirector('info.txt', 'warn.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = pet.AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from this file
    # (TODO: a link to raw data formats document to be given here)
    print('raw data: %s' % raw_data_file)
    acq_data = pet.AcquisitionData(raw_data_file)

    if ai_file is not None:
        anatomical_image = pet.ImageData()
        anatomical_image.read_from_file(ai_file)
        image = anatomical_image.get_uniform_copy()
        kappa = anatomical_image.get_uniform_copy()
        prior = pet.PLSPrior()
        prior.set_anatomical_image(anatomical_image)
        prior.set_kappa(kappa)
        prior.set_only_2D(True)
        print('using PLS prior...')
        aimage = prior.get_anatomical_image()
        aimage.show()
        tf = prior.get_only_2D()
        print(tf)
    else:
        prior = pet.QuadraticPrior()
        print('using Quadratic prior...')
        # create initial image estimate of dimensions and voxel sizes
        # compatible with the scanner geometry (included in the AcquisitionData
        # object ad) and initialize each voxel to 1.0
        image = acq_data.create_uniform_image(1.0)

    prior.set_up(image)
    prior.set_penalisation_factor(float(pen_factor))

    acq_model.set_up(acq_data, image)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = pet.make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)
    obj_fun.set_prior(prior)

    # select Ordered Subsets Maximum A-Posteriori One Step Late as the
    # reconstruction algorithm 
    recon = pet.OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)
    recon.set_input(acq_data)

    # set up the reconstructor based on the initial image
    # (checks the validity of parameters, sets up objective function
    # and other objects involved in the reconstruction, which involves
    # computing/reading sensitivity image etc etc.)
    print('setting up, please wait...')
    recon.set_up(image)

    # reconstruct from the initial image estimate
    # (check the OSEM demo to learn how to display results during sub-iterations)
    print('reconstructing...')
    recon.reconstruct(image)
    if show_plot:
        image.show(title = 'Reconstructed images')


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
