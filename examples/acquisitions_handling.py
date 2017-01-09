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

    # Get number of acquisitions:
    # the raw k-space data is a list of different 1D 
    # acquisitions (readouts) of different data type (e.g. noise correlation 
    # data, navigator data, image data,...).
    # The number of all aquisitions is
    na = input_data.number_of_acquisitions()
    # The number of image data acquisitions is
    ni = input_data.number_of_acquisitions('image')
    print('acquisitions: total %d, image data %d' % (na, ni))

    # sort data acquisition
    input_data.sort()

    where = range(254, 258)
    # inspect some acquisitions flags
    flags = input_data.get_info('flags')
    if flags[0] & IMAGE_DATA_MASK:
        print('first acquisition is image data')
    else:
        # should see this if input data file is test_2D_2x.h5
        print('first acquisition is not image data')
        
    # display flags for readout number 254 to 257    
    print('Flags'),
    print(flags[where])
    
    # inspect some kspace_encode_step_1 counters
    es1 = input_data.get_info('encode_step_1')
    print('Ky/PE - encoding'),
    print(es1[where])
    
    # inspect some slice counters
    slices = input_data.get_info('slice')
    print('Slices'),
    print(slices[where])
    
    # inspect some repetition counters
    print('Repetitions'),
    reps = input_data.get_info('repetition')
    print(reps[where])

    # copy raw data into python array and determine its size
    # in the case of the provided dataset 'simulated_MR_2D_cartesian.h5' the 
    # size is 2x256 phase encoding, 8 receiver coils and points 512 readout 
    # points (frequency encoding dimension)
    input_array = input_data.as_array()
    input_shape = input_array.shape
    print('input data dimensions: %dx%dx%d' % input_shape)

    # pre-process acquisitions
    print('---\n pre-processing acquisitions...')
    processed_data = PreprocessAcquisitions(input_data)

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
