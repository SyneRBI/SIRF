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

    # pre-process acquisitions
    print('pre-processing acquisitions...')
    preprocessed_data = input_data.process(['NoiseAdjustGadget', \
         'AsymmetricEchoGadget', 'RemoveROOversamplingGadget'])

    grappa_par = ''
    if args.no_gfactors:
        # gfactors are not needed
        grappa_par = '(send_out_gfactor=false)'

    # perform reconstruction
    nd = 'n_dimension=contrast,'
    sd = 's_dimension=average,'
    ss = 'split_slices=false,'
    ig = 'ignore_segment=true,'
    vb = 'verbose=true)'
    recon = ImagesReconstructor(['AcquisitionAccumulateTriggerGadget', \
         'BucketToBufferGadget(' + nd + sd + ss + ig + vb, \
         'PrepRefGadget', \
         'CartesianGrappaGadget' + grappa_par, \
         'FOVAdjustmentGadget', 'ScalingGadget', 'ImageArraySplitGadget'])
    recon.set_input(preprocessed_data)
    print('reconstructing...')
    recon.process()
    complex_output = recon.get_output()

    # post-process reconstructed images
    print('processing images...')
    complex_output.conversion_to_real(ISMRMRD_IMTYPE_MAGNITUDE)
    output = complex_output.process\
             (['ComplexToFloatGadget', 'FloatToShortGadget'])

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
