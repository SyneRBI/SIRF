'''OSMAPOSL reconstruction demo with user-controlled iterations

Usage:
  osmaposl_reconstruction [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -s <subs>, --subs=<subs>    number of subsets [default: 12]
  -i <iter>, --iter=<iter>    number of iterations [default: 2]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
  -p <path>, --path=<path>    sub-path to engine module [default: /xSTIR/pSTIR]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys
sys.path.append(os.environ.get('SRC_PATH') + args['--path'])
exec('from p' + args['--engine'] + ' import *')

# a simplistic example of user's involvement in the reconstruction
def my_image_data_processor(image_array, im_num):
    # plot the current estimate of the image at z = 20
    pylab.figure(im_num)
    pylab.title('image estimate %d' % im_num)
    pylab.imshow(image_array[20,:,:])
    # image is not modified in this simplistic example - but might have been
    return image_array

def main():

    # direct all engine's information and warnings printing to files
    printer = Printer('info.txt', 'warn.txt')

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())

    # PET acquisition data to be read from this file
    # (TODO: a link to raw data formats document to be given here)
    raw_data_file = args['--file']
    print('raw data: %s' % raw_data_file)
    ad = AcquisitionData(raw_data_file)

    # create initial image estimate of dimensions and voxel sizes
    # compatible with the scanner geometry (included in the AcquisitionData
    # object ad) and initialize each voxel to 1.0
    image = ad.create_empty_image(1.0)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood with linear model for mean
    # (TODO: find a good descriptive name for this object)
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    # select Ordered Subsets Maximum A-Posteriori One Step Late
    # as the reconstruction algorithm (since we are not using a penalty,
    # or prior, in this example, we will actually run OSEM)
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(int(args['--subs']))

    # set up the reconstructor based on a sample image
    # (checks the validity of parameters, sets up objective function
    # and other objects involved in the reconstruction, which involves
    # computing/reading sensitivity image etc etc.)
    print('setting up, please wait...')
    recon.set_up(image)

    # set the initial image estimate
    recon.set_current_estimate(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    for iteration in range(int(args['--iter'])):
        print('\n------------- iteration %d' % iteration)
        # perform one OSMAPOSL iteration
        recon.update_current_estimate()
        # copy current image estimate into python array to inspect/process
        image_array = recon.get_current_estimate().as_array()
        # apply user defined image data processor/visualizer
        my_image_array = my_image_data_processor(image_array, iteration + 1)
        # fill the current image estimate with new data
        recon.get_current_estimate().fill(my_image_array)
    pylab.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
