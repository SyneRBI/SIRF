'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by creating and running multiple gadget chains
of 3 types:
- acquisition processing chain
- reconstruction chain
- image processing chain
and how to visualise or modify data in between these chains.
'''

import argparse
import os
import sys
import matplotlib.pyplot as plt

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by creating and running multiple gadget chains
of 3 types:
- acquisition processing chain
- reconstruction chain
- image processing chain
and how to visualise or modify data in between these chains.
''')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

def gaussian(x, mu, sig):
    return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sig, 2.)))
    
def main():

    # Acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)
    # Get number of acquisitions
    # TODO: rename to number_of_acquisitions()
    na = input_data.number()
    
    # Get size of current k-space slices data as tuple
    # (number of samples, number of acquisitions per slice, number of coils)
    # TODO: reverse the tuple elements' order
    # TODO: replace with dimensions() method that returns the shape of as_array()
    kdim = input_data.slice_dimensions()
    # This way of printing works for both Python 2.* and Python 3.*
    print('Size of k-space slice reduced from %dx%dx%d' % kdim)
##    print "Size of k-space reduced from", kdim,
    
    # Pre-process acquisitions
    # Create an object which removes the readout oversampling from the acquired 
    # k-space data
    acq_proc = AcquisitionsProcessor(['RemoveROOversamplingGadget'])
    preprocessed_data = acq_proc.process(input_data)
    
    # Get size of k-space slices after removal of oversampling
    kdim = preprocessed_data.slice_dimensions()
    print('to %dx%dx%d' % kdim)
##    print "to", kdim, "."
    
    # Create simple Gaussian weighting function and apply it along the
    # readout direction onto the k-space data
    print('Apply Gaussian weighting function along readout')
##    print "Apply Gaussian weighting function along readout"
    gauss_weight = gaussian(numpy.array([numpy.linspace(-kdim[0]/2, kdim[0]/2, kdim[0])]),0,20)

    # This looks like the correct way of applying Gaussian weights
    gauss_weight = numpy.tile(gauss_weight, (na, 1))
##    gauss_weight = numpy.tile(gauss_weight.T, (1,kdim[1]))
    # TODO: make as_array() return ndarray of shape
    # (number of coils, number of slices, acquis. per slice, number of samples)
    # or even more dimensions
    dat_array = preprocessed_data.as_array()
##    dat_array = preprocessed_data.as_array()
    for c in range(kdim[2]):
        dat_array[:,c,:] = numpy.multiply(dat_array[:,c,:], gauss_weight)
##        dat_array[:,:,c] = numpy.multiply(dat_array[:,:,c], gauss_weight)

    # TODO: preprocessed_data.fill(dat_array)
    # Difficulty: need to make room for non-standard situations -
    # acquisitions of different type and shape etc.

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
