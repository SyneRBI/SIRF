'''
Upper-level interface demo that illustrates how MR data can be interfaced 
from python.

Usage:
  acquisitions_handling.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -r <rnge>, --range=<rnge>   range of acquisitions to examine as string '(a,b)'
                              [default: (254, 258)]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from ast import literal_eval
import os
import sys

# import engine module
exec('from p' + args['--engine'] + ' import *')

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    input_file = existing_filepath(data_path, args['--file'])

    # acquisitions will be read from an HDF file input_file
    input_data = AcquisitionData(input_file)

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

    # retrieve the range of acquisitions to examine
    t = literal_eval(args['--range'])
    if t[0] >= t[1] or t[1] >= na:
        raise error('Wrong acquisitions range')
    where = range(t[0], t[1])

    # retrieve acquisitions flags
    flags = input_data.get_info('flags')

    # inspect the first acquisition flag
    if flags[0] & IMAGE_DATA_MASK:
        print('first acquisition is image data')
    else:
        # should see this if input data file is test_2D_2x.h5
        print('first acquisition is not image data')
        
    # display flags
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
    processed_data = preprocess_acquisitions(input_data)

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
