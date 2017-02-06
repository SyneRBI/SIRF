'''
Lower-level interface demo that illustrates creating and running a chain
of gadgets - shortest version.

Usage:
  fully_sampled_recon_single_chain.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of $SRC_PATH/SIRF
  -o <file>, --output=<file>  images output file
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys

# locate the input data file
data_path = args['--path']
if data_path is None:
    SRC_PATH = os.environ.get('SRC_PATH')
    if SRC_PATH is None:
        print('Path to raw data files not set, please use -p <path> or --path=<path> to set it')
        sys.exit()
    data_path =  SRC_PATH + '/SIRF/data/examples/MR'
input_file = data_path + '/' + args['--file']
if not os.path.isfile(input_file):
    print('file %s not found' % input_file)

# import engine module
from pGadgetron import *

def main():

    # acquisitions will be read from an HDF file
    input_data = AcquisitionData(input_file)
    # create reconstruction object
    recon = ImagesReconstructor(['RemoveROOversamplingGadget', \
        'SimpleReconGadgetSet'])
    # reconstruct images
    images = recon.reconstruct(input_data)
    # show reconstructed images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
