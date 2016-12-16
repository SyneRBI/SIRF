'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by directly creating and running a chain of 
gadgets.
'''

import argparse
import os
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by directly creating and running a chain of 
gadgets.
''')
parser.add_argument('-o', '--output', default = None, help = 'output file name')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)
    
    # create reconstruction object
    # Rather than using a predefined image reconstruction object, here a new 
    # image reconstruction object is created by concatinating multiple gadgets 
    # (for more information on Gadgetron and its gadgets please see: 
    # https://github.com/gadgetron/.).
    # Parameters for individual gadgets can be defined either during the 
    # creation of the reconstruction object:
    #   e.g. AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)
    #   
    # or by using set_gadget_property()
    # The gadgets will be concatinated and will be executed as soon as 
    # process() is called
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

    if args.output is not None:
        # write images to a new group in args.output
        # named after the current date and time
        print('writing to %s' % args.output)
        time_str = time.asctime()
        images.write(args.output, time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
