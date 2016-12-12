'''
Lower-level interface demo that illustrates creating and running gadget chains
of 3 types:
- acquisition processing chain
- reconstruction chain
- image processing chain
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
Lower-level interface demo that illustrates creating and running gadget chains
of 3 types:
(i) acquisition processing chain
(ii) reconstruction chain
(iii) image processing chain
''')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)

    # pre-process acquisitions
    acq_proc = AcquisitionsProcessor(['RemoveROOversamplingGadget'])
    print('processing acquisitions...')
    preprocessed_data = acq_proc.process(input_data)

    # create reconstruction object
    recon = ImagesReconstructor\
        (['SimpleReconGadgetSet(trigger_dimension=repetition,split_slices=true)'])
    # connect to preprocessed data
    recon.set_input(preprocessed_data)
    # perform reconstruction
    print('reconstructing...')
    recon.process()
    # get reconstructed images
    complex_images = recon.get_output()

    # post-process reconstructed images
    img_proc = ImagesProcessor(['ExtractGadget'])
    print('processing images...')
    images = img_proc.process(complex_images)

    # show obtained images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
