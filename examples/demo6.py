'''
Lower-level demo, 3-chain GRAPPA reconstruction of undersampled data.
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
Lower-level demo, 3-chain GRAPPA reconstruction of undersampled data.
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
parser.add_argument('-o', '--output', default = None, help = 'output file name')
parser.add_argument\
('--no_gfactors', help = 'no gfactors to be computed', action = 'store_true')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)

    # define gadgets
    gadget11 = Gadget('NoiseAdjustGadget')
    gadget12 = Gadget('AsymmetricEchoGadget')
    gadget13 = Gadget('RemoveROOversamplingGadget')
    gadget21 = Gadget('AcquisitionAccumulateTriggerGadget')
    gadget22 = Gadget('BucketToBufferGadget')
    gadget23 = Gadget('PrepRefGadget')
    gadget24 = Gadget('CartesianGrappaGadget')
    gadget25 = Gadget('FOVAdjustmentGadget')
    gadget26 = Gadget('ScalingGadget')
    gadget27 = Gadget('ImageArraySplitGadget')
    gadget31 = Gadget('ComplexToFloatGadget')
    gadget32 = Gadget('FloatToShortGadget')

    if args.no_gfactors:
        # gfactors are not needed
        gadget24.set_property('send_out_gfactor', 'false')

    # pre-process acquisitions
    acq_proc = AcquisitionsProcessor()
    acq_proc.add_gadget('g1', gadget11)
    acq_proc.add_gadget('g2', gadget12)
    acq_proc.add_gadget('g3', gadget13)
    print('pre-processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)

    # create reconstruction object
    recon = ImagesReconstructor()
    recon.add_gadget('g1', gadget21)
    recon.add_gadget('g2', gadget22)
    recon.add_gadget('g3', gadget23)
    recon.add_gadget('g4', gadget24)
    recon.add_gadget('g5', gadget25)
    recon.add_gadget('g6', gadget26)
    recon.add_gadget('g7', gadget27)
    # connect to input data
    recon.set_input(preprocessed_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    complex_output = recon.get_output()

    # post-process reconstructed images
    img_proc = ImagesProcessor()
    img_proc.add_gadget('g1', gadget31)
    img_proc.add_gadget('g2', gadget32)
    complex_output.conversion_to_real(1)
    print('processing images...')
    output = img_proc.process(complex_output)

    # show images
    output.show()

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
    print ('Gadgetron exception occured:\n', err.value)
