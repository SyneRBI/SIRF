'''
Demonstrates use of the EPSRC-funded CCP-PETMR code (SIRF). 
See function grappa_detail.m for an example showing more of the 
workings and functionality of the SIRF code.

Pre-requisites:
 1) This MATLAB code needs to be able to access a listening gadgetron.
    On the Virtual Machine, gadgetron is installed and the user just needs
    to type 'gadgetron' in a terminal window.
    On standalone systems, the user will need to have installed ISMRMRD
    and gadgetron code.

 2) An input data file from a GRAPPA MRI acquisition in the ISMRMRD format.
    Example GRAPPA datasets:
    a) 'meas_MID00108_FID57249_test_2D_2x.dat' is 
       available from https://www.ccppetmr.ac.uk/downloads
       This is in the manufacturer's raw data format and needs to be
       converted to ISMRMRD format using 'siemens_to_ismrmrd'.
       This executable is installed on the Virtual Machine.

    b) A simulated ISMRMRD h5 file is available as default

Usage:
  grappa_basic.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
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
    input_file = existing_filepath(data_path, args['--file'])
    
    # Initially we create a container that points to the h5 file. Data is
    # not read from file until the gadgetron is called using
    # the 'process' method.
    
    # Create an Acquisition Container of type pGadgetron.AcquisitionData
    print('---\n reading in file %s...' % input_file)
    input_Cont = AcquisitionData(input_file)
    
    
    # Pre-process this input data. (Currently this is a MATLAB script that just
    # sets up a 3 chain gadget. In the future it will be independent of the MR
    # recon engine.)
    print('---\n pre-processing acquisitions...')
    preprocessed_AcCont = preprocess_acquisitions(input_Cont)
    
    
    # Perform reconstruction of the preprocessed data.
    # 1. set the reconstruction to be for Cartesian GRAPPA data.
    recon = GenericCartesianGRAPPAReconstruction();
    
    # 2. set the reconstruction input to be the data we just preprocessed.
    recon.set_input(preprocessed_AcCont);
    
    # 3. run (i.e. 'process') the reconstruction.
    print('---\n reconstructing...\n');
    recon.process();
    

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