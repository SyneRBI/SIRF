'''OSMAPOSL reconstruction demo with user-controlled iterations

Usage:
  osmaposl_reconstruction [--help | options]

Options:
  -e=<e>, --engine=<e>  reconstruction engine [default: Stir]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys
sys.path.append(os.environ.get('SRC_PATH') + '/xSTIR/pSTIR')
exec('from p' + args['--engine'] + ' import *')

# a simplistic example of user's involvement in the reconstruction
def my_image_data_processor(image_array, im_num):
    # plot the current estimate of the image
    pylab.figure(im_num)
    pylab.title('image estimate %d' % im_num)
    pylab.imshow(image_array[20,:,:])
    print('close Figure %d window to continue' % im_num)
    pylab.show()
    # image is not modified in this simplistic example - but might have been
    return image_array

def main():

    # direct all information printing and warnings to files
    info_printer = printerTo('info.txt', INFO_CHANNEL)
    warning_printer = printerTo('warn.txt', WARNING_CHANNEL)

    # select acquisition model that implements the geometric
    # forward projection by a ray tracing matrix multiplication
    am = AcquisitionModelUsingMatrix(RayTracingMatrix())

    # PET acquisition data to be read from this file
    ad = AcquisitionData('my_forward_projection.hs')

    # create initial image estimate of dimensions and voxel sizes
    # compatible with acquisition data and initialize each voxel to 1.0
    image = ad.create_empty_image(1.0)

    # define objective function to be maximized as
    # Poisson logarithmic likelihood with linear model for mean
    obj_fun = PoissonLogLh_LinModMean_AcqMod()
    obj_fun.set_acquisition_model(am)
    obj_fun.set_acquisition_data(ad)

    # select Ordered Subsets Maximum A-Posteriori One Step Late
    # as the reconstruction algorithm (since we are not using a penalty
    # (or prior) in this example, we will actually run OSEM)
    recon = OSMAPOSLReconstruction()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(12)

    # set up the reconstructor
    print('setting up, please wait...')
    recon.set_up(image)

    # in order to see the reconstructed image evolution
    # take over the control of the iterative process
    # rather than allow recon.reconstruct to do all job at once
    recon.set_current_estimate(image)
    for iteration in range(2):
        print('\n------------- iteration %d' % iteration)
        # perform an iteration
        recon.update_current_estimate()
        # copy current image estimate into python array to inspect/process
        image_array = recon.get_current_estimate().as_array()
        # apply your image data processor/visualizer
        my_image_array = my_image_data_processor(image_array, iteration + 1)
        # fill the current image estimate with new data
        recon.get_current_estimate().fill(my_image_array)

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
except error as err:
    # display error information
    print('STIR exception occured: %s' % err.value)
