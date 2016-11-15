'''
Upper-level demo, reconstruction of arbitrarily sampled data.
'''

import argparse
import os
import pylab
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Upper-level demo, reconstruction of arbitrarily sampled data.
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoAdjustROGadget', \
         'RemoveROOversamplingGadget']
    print('---\n pre-processing acquisitions...')
    preprocessed_data = input_data.process(prep_gadgets)

    # perform reconstruction
    if not input_data.is_undersampled():
        print('---\n reconstructing fully sampled data...')
        recon = MR_BasicReconstruction()
    else:
        print('---\n reconstructing undersampled data using GRAPPA...')
        recon = MR_BasicGRAPPAReconstruction()
        recon.compute_gfactors(False)
    complex_images = recon.reconstruct(preprocessed_data)
    images = complex_images.real()

    # show reconstructed images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
