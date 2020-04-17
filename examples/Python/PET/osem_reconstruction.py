'''OSEM reconstruction demo.
We actually use the OSMAPOSL reconstructor in this demo. This reconstructor
implements an Ordered Subsets (OS) version of the One Step Late algorithm (OSL)
from Green et al for Maximum a Posteriori (MAP) maximisation. Here we use it
for Maximum Likelihood (ML) in which case it is equivalent to OSEM.

Usage:
  osem_reconstruction [--help | options]

Options:
  -f <file>, --file=<file>     raw data file [default: my_forward_projection.hs]
  -a <file>, --anim=<file>     anatomical image file
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
  -s <subs>, --subs=<subs>     number of subsets [default: 12]
  -i <siter>, --subiter=<siter>  number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
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

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

# process command-line options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
if args['--anim'] is not None:
    ai_file = existing_filepath(data_path, args['--anim'])
else:
    ai_file = None

# Define a function that does something with an image. This function
# provides a simplistic example of user's involvement in the reconstruction
def image_data_processor(image_array, im_num):
    """ Process/display an image"""
    # display the current estimate of the image at z = 20
    pylab.figure(im_num)
    pylab.title('image estimate %d' % im_num)
    pylab.imshow(image_array[20,:,:])
    print('You may need to close Figure %d window to continue' % im_num)
    # image is not modified in this simplistic example - but might have been
    return image_array

def main():
 
    # direct all engine's information and warnings printing to files
    msg_red = MessageRedirector('info.txt', 'warn.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    acq_model = AcquisitionModelUsingRayTracingMatrix()

    # PET acquisition data to be read from this file
    # (TODO: a link to raw data formats document to be given here)
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # create initial image estimate of dimensions and voxel sizes
    # compatible with the scanner geometry (included in the AcquisitionData
    # object ad) and initialize each voxel to 1.0
    image = acq_data.create_uniform_image(1.0)

    acq_model.set_up(acq_data, image)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = make_Poisson_loglikelihood(acq_data)
    obj_fun.set_acquisition_model(acq_model)

    # select Ordered Subsets Maximum A-Posteriori One Step Late as the
    # reconstruction algorithm (since we are not using a penalty, or prior, in
    # this example, we actually run OSEM);
    # this algorithm does not converge to the maximum of the objective function
    # but is used in practice to speed-up calculations
    recon = OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_input(acq_data)

    # set up the reconstructor based on a sample image
    # (checks the validity of parameters, sets up objective function
    # and other objects involved in the reconstruction, which involves
    # computing/reading sensitivity image etc etc.)
    print('setting up, please wait...')
    recon.set_up(image)

    # set the initial image estimate
    recon.set_current_estimate(image)

    # in order to see the reconstructed image evolution
    # open up the user's access to the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for subiteration in range(num_subiterations):
        print('\n------------- sub-iteration %d' % subiteration)
        # perform one OSMAPOSL sub-iteration
        recon.update_current_estimate()
        # copy current image estimate into python array to inspect/process
        image_array = recon.get_current_estimate().as_array()
        # apply user defined image data processor/visualizer
        processed_image_array = image_data_processor(image_array, subiteration + 1)
        # fill the current image estimate with new data
        image.fill(processed_image_array)
        recon.set_current_estimate(image)
    pylab.show()

    # forward projection of the reconstructed image simulates the
    # acquisition of data by the scanner
    print('projecting...')
    simulated_data = acq_model.forward(image)
    # compute the reconstruction residual
    diff = simulated_data * (acq_data.norm()/simulated_data.norm()) - acq_data
    print('relative residual norm: %e' % (diff.norm()/acq_data.norm()))

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
