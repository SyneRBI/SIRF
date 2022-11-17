'''
Medium-level interface demo that illustrates Cartesian MR image
reconstruction using Gadgetron by creating and running several gadget chains
of different kind:
- acquisition processing chain
- reconstruction chain
- image processing chain
and how to visualise or modify data in between these chains.

Usage:
  recon_by_several_chains.py [--help | options]

Options:
  -f <file>, --file=<file>     raw data file
                               [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/MR
                               subfolder of SIRF root folder
  -s <sigma>, --sigma=<sigma>  gaussian sigma [default: 20]
  -o <outp>, --output=<path>   output file name [default: images.h5]
  --non-interactive            do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF).
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC.
## Copyright 2015 - 2017 University College London.
## Copyright 2015 - 2017 Physikalisch-Technische Bundesanstalt.
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

import numpy

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module objects
from pGadgetron import examples_data_path
from pGadgetron import existing_filepath
from pGadgetron import AcquisitionData
from pGadgetron import AcquisitionDataProcessor
from pGadgetron import Reconstructor
from pGadgetron import ImageDataProcessor
from pGadgetron import error

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
sigma = float(args['--sigma'])
show_plot = not args['--non-interactive']


def gaussian(x, mu, sigma):
    return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sigma, 2.)))


def main():

    # Acquisitions will be read from this HDF file
    input_file = existing_filepath(data_path, data_file)
    acq_data = AcquisitionData(input_file)
    undersampled = acq_data.is_undersampled()

    # Create and run Gadgetron chain that pre-processes acquisition data
    # (removes the readout oversampling from the acquired k-space data etc.)
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoAdjustROGadget', \
                                'RemoveROOversamplingGadget' ]
    acq_proc = AcquisitionDataProcessor(prep_gadgets)
    acq_proc.set_input(acq_data)
    acq_proc.process()
    preprocessed_data = acq_proc.get_output()
    # shortcut for the above 3 lines
##    preprocessed_data = acq_proc.process(acq_data)

    # Compare sizes of k-space data as tuples
    # (number of acquisitions, number of coils, number of samples)
    # before and after removal of oversampling
    k_space_dimensions = preprocessed_data.dimensions()
    print('Size of k-space slice reduced from %dx%dx%d' % acq_data.dimensions())
    print('to %dx%dx%d' % k_space_dimensions)
    
    # Create simple Gaussian weighting function and apply it along the
    # readout direction onto the k-space data
    print('Apply Gaussian weighting function along readout')
    kdim = k_space_dimensions
    gauss_weight = gaussian\
        (numpy.array([numpy.linspace(-kdim[2]/2, kdim[2]/2, kdim[2])]),0,sigma)
    gauss_weight = numpy.tile(gauss_weight, (kdim[0], 1))
    preprocessed_array = preprocessed_data.as_array()
    for c in range(kdim[1]):
        preprocessed_array[:,c,:] = \
            numpy.multiply(preprocessed_array[:,c,:], gauss_weight)
    preprocessed_data.fill(preprocessed_array)

    # create reconstruction chain
    recon_gadgets = [
        'GenericReconCartesianReferencePrepGadget',
        'GRAPPA:GenericReconCartesianGrappaGadget',
        'GenericReconFieldOfViewAdjustmentGadget',
        'GenericReconImageArrayScalingGadget'] if undersampled else ['SimpleReconGadget']
    recon = Reconstructor \
        (['AcquisitionAccumulateTriggerGadget', 'BucketToBufferGadget']
         + recon_gadgets + ['ImageArraySplitGadget'])
    
    # provide pre-processed k-space data
    recon.set_input(preprocessed_data)
    
    # perform reconstruction
    recon.process()

    # retrieve the reconstructed complex images
    if undersampled:
        reconstructed_data = recon.get_output('image')
    else:
        reconstructed_data = recon.get_output()

    # convert the images to real and scale using images processing chain
    img_proc = ImageDataProcessor(['ExtractGadget', 'AutoScaleGadget'])
    img_proc.set_input(reconstructed_data)
    img_proc.process()
    image_data = img_proc.get_output()
    print(image_data.dimensions())

    # write the images to file;
    # if the file name extension is .dcm, the writing is in DICOM format
    # and performed by yet another image processing chain
    image_data.write(args['--output'])

    # show obtained images
    if show_plot:
        image_data.show(title = 'Reconstructed image data (magnitude)')


try:
    main()
    print('\n=== done with %s' % __file__)
except error as err:
    # display error information
    print('??? %s' % err.value)
    #exit(1)
