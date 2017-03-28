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
                              subfolder of SIRF root folder
  -s=<sigma>, --sigma=<sigma>  gaussian sigma [default: 20]
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

sigma = float(args['--sigma'])

# import engine module
from pGadgetron import *

def gaussian(x, mu, sigma):
    return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sigma, 2.)))
    
def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    # Acquisitions will be read from this HDF file
    input_file = existing_filepath(data_path, args['--file'])
    input_data = AcquisitionData(input_file)

    # Get size of current k-space data as tuple
    # (number of acquisitions, number of coils, number of samples)
    kdim = input_data.dimensions()
    # This way of printing works for both Python 2.* and Python 3.*
    print('Size of k-space slice reduced from %dx%dx%d' % kdim)
    
    # Pre-process acquisitions
    # Create an object which removes the readout oversampling from the acquired 
    # k-space data
    acq_proc = AcquisitionDataProcessor(['RemoveROOversamplingGadget'])
    preprocessed_data = acq_proc.process(input_data)
    
    # Get size of k-space data after removal of oversampling
    kdim = preprocessed_data.dimensions()
    print('to %dx%dx%d' % kdim)
    
    # Create simple Gaussian weighting function and apply it along the
    # readout direction onto the k-space data
    print('Apply Gaussian weighting function along readout')
    gauss_weight = gaussian\
        (numpy.array([numpy.linspace(-kdim[2]/2, kdim[2]/2, kdim[2])]),0,sigma)
    gauss_weight = numpy.tile(gauss_weight, (kdim[0], 1))
    data_array = preprocessed_data.as_array()
    for c in range(kdim[1]):
        data_array[:,c,:] = numpy.multiply(data_array[:,c,:], gauss_weight)

    preprocessed_data.fill(data_array)

    # create reconstruction object
    recon = Reconstructor\
        (['AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)', \
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
    img_proc = ImageDataProcessor(['ExtractGadget'])
    images = img_proc.process(complex_images)

    # show obtained images
    image_array = images.as_array()
    title = 'Reconstructed images (magnitude)'
    show_3D_array(abs(image_array), suptitle = title, label = 'slice')

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
