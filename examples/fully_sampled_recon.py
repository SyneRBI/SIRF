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
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

exec('from ' + args.engine + ' import *')

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)

    # pre-process acquisition data
    print('---\n pre-processing acquisitions...')
    processed_data = MR_remove_x_oversampling(input_data)

    # perform reconstruction
    recon = SimpleReconstruction()
    recon.set_input(processed_data)
    print('---\n reconstructing...')
    recon.process()
    images = recon.get_output()

    # show obtained images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
