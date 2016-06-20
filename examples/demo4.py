'''
Upper-level interface demo that illustrates pre-processing of acquisitions,
reconstructing images and post-processing them.
'''

import argparse
import os
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Upper-level interface demo that illustrates pre-processing of acquisitions,
reconstructing images and post-processing them.
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)

    # pre-process acquisition data
    print('processing acquisitions...')
    processed_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    print('reconstructing...')
    recon.process()
    complex_images = recon.get_output()

    # post-process reconstructed images
    print('processing images...')
    images = MR_extract_real_images(complex_images)

    # show obtained images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
