'''
Example of radial phase encoding (RPE) reconstruction.

Upper-level demo that illustrates the computation of how to use a non-cartesian 
radial phase-encoding acquisition model to reconstruct data. The computed
density compensation function simply accounts for multiply acquired points.

Usage:
  rpe_recon.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: 3D_RPE_Lowres.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
  -n <bool>, --non-cartesian  run recon iff non-cartesian code was compiled 
                              [default: false]
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
## (http://www.ccpsynerbi.ac.uk/).
##
## Licensed under the Apache License, Version 2.0 (the "License");
##   you may not use this file except in compliance with the License.
##   You may obtain a copy of the License at
##       http://www.apache.org/licenses/LICENSE-2.0
##   Unless required by applicable law or agreed to in writing, software
##   distributed under the License is distributed on an "AS IS" BASIS,
##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
##   See the License for the specific language governing permissions and
##   limitations under the License.

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR') + '/zenodo/'
output_file = args['--output']
show_plot = not args['--non-interactive']
run_recon = args['--non-cartesian']

import numpy as np
    
def calc_unit_dcf(acq_data):

    traj = np.transpose(get_grpe_trajectory(acq_data))

    traj, inverse, counts = np.unique(traj, return_inverse=True, return_counts=True, axis=1)
    
    dcf = ( 1.0 / counts)[inverse]
    max_traj_rad = np.max(np.linalg.norm(traj, axis=0))
    dcf_norm =  np.sum(dcf) / (max_traj_rad**2 * np.pi)
    dcf = dcf / dcf_norm
  
    return dcf


def main():

    # locate the k-space raw data file
    input_file = existing_filepath(data_path, data_file)
    
    # acquisition data will be read from an HDF file input_file
    AcquisitionData.set_storage_scheme('memory')
    acq_data = AcquisitionData(input_file)
    
    print('---\n acquisition data norm: %e' % acq_data.norm())

    # pre-process acquisition data
    print('---\n pre-processing acquisition data...')
    processed_data  = preprocess_acquisition_data(acq_data)
    
    # sort processed acquisition data;
    print('---\n sorting acquisition data...')
    processed_data.sort()
    
    #set the trajectory and compute the dcf
    print('---\n setting the trajectory...')
    processed_data = set_grpe_trajectory(processed_data)
    
    print('---\n computing density weights...')
    dcf = calc_unit_dcf(processed_data)
    processed_data = set_densitycompensation_as_userfloat(processed_data, dcf)

    print("Am i running the rest of the code? : " + str(run_recon))
    
    if run_recon is True:
        print('---\n computing coil sensitivity maps...')
        csms = CoilSensitivityData()
        csms.smoothness = 10
        
        csms.calculate(processed_data)
        
        # create acquisition model based on the acquisition parameters
        print('---\n Setting up Acquisition Model...')
    
        acq_model = AcquisitionModel()
        acq_model.set_up(processed_data, csms.copy())
        acq_model.set_coil_sensitivity_maps(csms)
    
        print('---\n Backward projection ...')
        recon_img = acq_model.backward(processed_data)
        
        if show_plot:
            recon_img.show(title = 'Reconstructed images (magnitude)')
    else:
        print('---\n Skipping non-cartesian code...')

    
        


try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
