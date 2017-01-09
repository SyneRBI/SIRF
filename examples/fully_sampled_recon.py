'''
Low-level interface demo that illustrates pre-processing of MR raw (k-
space) data, 2D image reconstruction using FFT and image display.
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
Low-level interface demo that illustrates pre-processing of MR raw (k-
space) data, 2D image reconstruction using FFT and image display.
''')
parser.add_argument('-e', '--engine', default = 'pGadgetron', help = 'engine')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

exec('from ' + args.engine + ' import *')

def main():

    # acquisitions will be read from an HDF file args.filename
    # MR raw data formats from different vendors can be transformed to 
    # HDF file format using siemens_to_ismrmrd, philips_to_ismrmrd or
    # bruker_to_ismrmrd on https://github.com/ismrmrd/.
    print('---\n reading in file %s...' % args.filename)
    input_data = AcquisitionData(args.filename)

    # pre-process acquired k-space data
    # Prior to image reconstruction several pre-processing steps such as 
    # assymetric echo compensation, noise decorelation for multi-coil data or 
    # removal of oversampling along frequency encoding (i.e. readout or kx)
    # direction. So far only the removal of readout oversampling and noise and
    # asymmetric echo adjusting is implemented
    print('---\n pre-processing acquisitions...')
    processed_data = PreprocessAcquisitions(input_data)

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
