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

#from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Upper-level interface demo that illustrates pre-processing of acquisitions,
reconstructing images and post-processing them.
''')
parser.add_argument('-e', '--engine', default = 'pGadgetron', help = 'engine')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

exec 'from ' + args.engine + ' import *'

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)

    # pre-process acquisition data
    print('processing acquisitions...')
    processed_data = input_data.process(['RemoveROOversamplingGadget'])

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    print('reconstructing...')
    recon.process()
    complex_images = recon.get_output()

    # extract real images from complex
    print('processing images...')
    images = complex_images.real()

    # show obtained images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
