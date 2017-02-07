'''
Low-level interface demo that illustrates pre-processing of MR raw (k-
space) data, 2D image reconstruction using FFT and image display.

Usage:
  fully_sampled_recon.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of $SRC_PATH/SIRF
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
import sys

# import engine module
exec('from p' + args['--engine'] + ' import *')

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    input_file = existing_file(data_path, args['--file'])

    # acquisitions will be read from an HDF file
    # MR raw data formats from different vendors can be transformed to 
    # HDF file format using siemens_to_ismrmrd, philips_to_ismrmrd or
    # bruker_to_ismrmrd on https://github.com/ismrmrd/.
    print('---\n reading in file %s...' % input_file)
    input_data = AcquisitionData(input_file)

    # pre-process acquired k-space data
    # Prior to image reconstruction several pre-processing steps such as 
    # assymetric echo compensation, noise decorelation for multi-coil data or 
    # removal of oversampling along frequency encoding (i.e. readout or kx)
    # direction. So far only the removal of readout oversampling and noise and
    # asymmetric echo adjusting is implemented
    print('---\n pre-processing acquisitions...')
    processed_data = preprocess_acquisitions(input_data)

    # setup reconstruction
    # Create a reconstruction object (in this case simple 2D Cartesian FFT) and
    # provide pre-processed k-space data as input
    recon = SimpleReconstruction()
    recon.set_input(processed_data)
    
    # perform reconstruction
    print('---\n reconstructing...')
    recon.process()
    
    # retrieve reconstruced images
    images = recon.get_output()

    # show reconstructed images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
