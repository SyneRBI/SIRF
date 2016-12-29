'''
Upper-level interface demo that illustrates how MR data can be interfaced 
from python.
'''

import argparse
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Upper-level interface demo that illustrates how MR data can be interfaced 
from python.
''')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)

    # sort data acquisition
    # prior to this step the raw k-space data is a list of different 1D 
    # acquisitions (readouts) of different data type (e.g. noise correlation 
    # data, navigator data, image data,...). Afterwards input_data contains 
    # only image data
    input_data.sort()

    # copy raw data into python array and determine its size
    # in the case of the provided dataset 'simulated_MR_2D_cartesian.h5' the 
    # size is 256 phase encoding, 8 receiver coils and points 512 readout 
    # points (frequency encoding dimension)
    input_array = input_data.as_array()
    input_shape = input_array.shape
    print('input data dimensions: %dx%dx%d' % input_shape)

    # remove oversampling along readout
    processed_data = input_data.process(['RemoveROOversamplingGadget'])

    # copy processed acquisitions into an array and determine its size
    # by removing the oversampling factor of 2 along the readout direction, the
    # number of readout samples was halfed
    processed_array = processed_data.as_array()
    processed_shape = processed_array.shape
    print('processed data dimensions: %dx%dx%d' % processed_shape)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
