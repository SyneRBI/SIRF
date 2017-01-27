'''
Lower-level interface demo that illustrates creating and running a chain
of gadgets - shortest version.

Usage:
  fully_sampled_recon_single_chain.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    sub-path to engine module
                              [default: /xGadgetron/pGadgetron]
  -o <file>, --output=<file>  images output file
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys
import time

sys.path.append(os.environ.get('SRC_PATH') + args['--path'])

from pGadgetron import *

def main():

    # acquisitions will be read from an HDF file
    input_data = AcquisitionData(args['--file'])
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
