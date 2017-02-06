'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by creating and running multiple gadget chains
of 3 types:
- acquisition processing chain
- reconstruction chain
- image processing chain
and how to visualise or modify data in between these chains.

Usage:
  fully_sampled_recon_three_chains.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of $SRC_PATH/SIRF
  -s=<sigma>, --sigma=<sigma>  gaussian sigma [default: 20]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys

# locate the input data file
data_path = args['--path']
if data_path is None:
    SRC_PATH = os.environ.get('SRC_PATH')
    if SRC_PATH is None:
        print('Path to raw data files not set, please use -p <path> or --path=<path> to set it')
        sys.exit()
    data_path =  SRC_PATH + '/SIRF/data/examples/MR'
input_file = data_path + '/' + args['--file']
if not os.path.isfile(input_file):
    print('file %s not found' % input_file)

sigma = float(args['--sigma'])

# import engine module
from pGadgetron import *

def gaussian(x, mu, sig):
    return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sig, 2.)))
    
def main():

    # Acquisitions will be read from an HDF file
    input_data = AcquisitionData(input_file)

    # Get size of current k-space data as tuple
    # (number of acquisitions, number of coils, number of samples)
    kdim = input_data.dimensions()
    # This way of printing works for both Python 2.* and Python 3.*
    print('Size of k-space slice reduced from %dx%dx%d' % kdim)
    
    # Pre-process acquisitions
    # Create an object which removes the readout oversampling from the acquired 
    # k-space data
    acq_proc = AcquisitionsProcessor(['RemoveROOversamplingGadget'])
    preprocessed_data = acq_proc.process(input_data)
    
    # Get size of k-space data after removal of oversampling
    kdim = preprocessed_data.dimensions()
    print('to %dx%dx%d' % kdim)
    
    # Create simple Gaussian weighting function and apply it along the
    # readout direction onto the k-space data
    print('Apply Gaussian weighting function along readout')
    gauss_weight = gaussian(numpy.array([numpy.linspace(-kdim[2]/2, kdim[2]/2, kdim[2])]),0,sigma)
    gauss_weight = numpy.tile(gauss_weight, (kdim[0], 1))
    data_array = preprocessed_data.as_array()
    for c in range(kdim[1]):
        data_array[:,c,:] = numpy.multiply(data_array[:,c,:], gauss_weight)

    preprocessed_data.fill(data_array)

    # create reconstruction object
    recon = ImagesReconstructor(['AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)', \
        'BucketToBufferGadget(split_slices=true, verbose=false)', 
        'SimpleReconGadget', 'ImageArraySplitGadget'])
    
    # provide pre-processed k-space data
    recon.set_input(preprocessed_data)
    
    # perform reconstruction
    recon.process()
    
    # retrieve reconstructed images
    complex_images = recon.get_output()

    # post-process reconstructed images
    # Rather than using the function get_output(), a new object based on the 
    # gadget "ExtractGadget" is specified to retrieve the image data
    img_proc = ImagesProcessor(['ExtractGadget'])
    images = img_proc.process(complex_images)

    # show obtained images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
