'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by directly creating and running a chain of 
gadgets.

Usage:
  fully_sampled_recon_single_chain.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  images output file
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import time

output_file = args['--output']

# import engine module
from pGadgetron import *

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    input_file = existing_filepath(data_path, args['--file'])

    # acquisitions will be read from an HDF file
    input_data = AcquisitionData(input_file)
    
    # create reconstruction object
    # Rather than using a predefined image reconstruction object, here a new 
    # image reconstruction object is created by concatinating multiple gadgets 
    # (for more information on Gadgetron and its gadgets please see: 
    # https://github.com/gadgetron/.).
    # Parameters for individual gadgets can be defined either during the 
    # creation of the reconstruction object:
    #   e.g. AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)
    # or by giving a gadget a label (cf. label ex: for the last gadget)
    # and using set_gadget_property(label, propery, value).
    # The gadgets will be concatenated and will be executed as soon as 
    # process() is called.
    recon = ImagesReconstructor(['RemoveROOversamplingGadget', \
        'AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)', \
        'BucketToBufferGadget(split_slices=true, verbose=false)', \
        'SimpleReconGadget', 'ImageArraySplitGadget', 'ex:ExtractGadget'])
        
    # ExtractGadget defines which type of image should be returned:
    # none      0
    # magnitude 1
    # real      2
    # imag      4
    # phase     8
    # max       16  
    # in this example '5' returns both magnitude and imag    
    recon.set_gadget_property('ex', 'extract_mask', 5) 
    
    # provide raw k-space data as input
    recon.set_input(input_data)
    
    # perform reconstruction
    recon.process()
    
    # retrieve reconstructed images
    images = recon.get_output()

    # show reconstructed images
    images.show()

    if output_file is not None:
        # write images to a new group in args.output
        # named after the current date and time
        time_str = time.asctime()
        print('writing to %s' % output_file)
        images.write(output_file, time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
