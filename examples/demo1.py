'''
Lower-level interface demo that illustrates creating and running a chain
of gadgets.
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
Lower-level interface demo that illustrates creating and running a chain
of gadgets.
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

def main():
    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)
    
    # define gadgets
    gadget1 = Gadget('RemoveROOversamplingGadget')
    gadget2 = Gadget('AcquisitionAccumulateTriggerGadget')
    gadget3 = Gadget('BucketToBufferGadget')
    gadget4 = Gadget('SimpleReconGadget')
    gadget5 = Gadget('ImageArraySplitGadget')
    gadget6 = Gadget('ExtractGadget')

    # set gadget properties
    gadget2.set_property('trigger_dimension', 'repetition')
    gadget3.set_property('split_slices', 'true')

    # create reconstruction object
    recon = ImagesReconstructor()

    # build gadgets chain
    recon.add_gadget('g1', gadget1)
    recon.add_gadget('g2', gadget2)
    recon.add_gadget('g3', gadget3)
    recon.add_gadget('g4', gadget4)
    recon.add_gadget('g5', gadget5)
    recon.add_gadget('g6', gadget6)

    # connect to input data
    recon.set_input(input_data)
    # perform reconstruction
    recon.process()
    
    # get reconstructed images
    images = recon.get_output()

    # show reconstructed images
    images.show()

    # write images to a new group in 'output1.h5'
    # named after the current date and time
    print('appending output1.h5...')
    time_str = time.asctime()
    images.write('output1.h5', time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
