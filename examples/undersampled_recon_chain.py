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
('--gfactors', help = 'gfactors to be computed', action = 'store_true')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)

    # pre-process acquisitions
    print('---\n pre-processing acquisitions...')
    preprocessed_data = input_data.process(['NoiseAdjustGadget', \
         'AsymmetricEchoAdjustROGadget', 'RemoveROOversamplingGadget'])

    # set up reconstruction chain
    recon = ImagesReconstructor([\
         'AcquisitionAccumulateTriggerGadget', 'BucketToBufferGadget' \
         '(N_dimension=contrast,S_dimension=average,split_slices=false)', \
         'GenericReconCartesianReferencePrepGadget', \
         'GRAPPA:GenericReconCartesianGrappaGadget', \
         'GenericReconFieldOfViewAdjustmentGadget', \
         'GenericReconImageArrayScalingGadget', 'ImageArraySplitGadget'])
    # change a property of the gadget labelled by 'GRAPPA'
    recon.set_gadget_property('GRAPPA', 'send_out_gfactor', args.gfactors)
    recon.set_input(preprocessed_data)
    # reconstruct
    print('---\n reconstructing...')
    recon.process()
    output = recon.get_output()

    # show images
    output.show()

    if args.output is not None:
        # write images to a new group in args.output
        # named after the current date and time
        print('writing to %s' % args.output)
        time_str = time.asctime()
        output.write(args.output, time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
