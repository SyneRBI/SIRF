'''
Upper-level demo that illustrates the computation of how to use a non-cartesian
radial phase-encoding acquisition model to reconstruct data. The computed 
density compensation function does not supply optimal weights.

Usage:
  rpe_recon.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: Lowres_RPE.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
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
    data_path = examples_data_path('MR')
output_file = args['--output']
show_plot = not args['--non-interactive']

import numpy as np
from scipy import stats


def calc_dcf(traj):
    # this is an educated guess of the DCF using scipy kde
    traj, inverse, counts = np.unique(traj, return_inverse=True, return_counts=True, axis=1)
        
    bandwith_method = 0.05
    kernel = stats.gaussian_kde(traj, bw_method=bandwith_method)
    dcf = ( 1 / kernel(traj).T / counts)[inverse] # down-weight double samples
        
    max_traj_rad = np.max(np.linalg.norm(traj, axis=0))
    dcf_norm =  np.sum(dcf) / (max_traj_rad**2 * np.pi)  
    dcf = dcf / dcf_norm
    
    if show_plot:
        kmin, kmax = -max_traj_rad , max_traj_rad 
        nsteps = 100j        
        KX,KY = np.mgrid[kmin:kmax:nsteps, kmin:kmax:nsteps]
    
        interpol_pos = np.vstack([KX.ravel(), KY.ravel()])    
        dcf_interpol = np.reshape(kernel(interpol_pos).T, KX.shape)
    
        import matplotlib.pyplot as plt
    
        fig, ax = plt.subplots()
        ax.imshow(np.rot90(dcf_interpol), cmap=plt.cm.gist_earth_r, extent=[kmin, kmax, kmin, kmax])
        ax.plot(traj[0,:], traj[1,:], 'k.', markersize=2)
        ax.set_xlim([kmin, kmax])
        ax.set_ylim([kmin, kmax])
        plt.show()
    
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
    processed_data = preprocess_acquisition_data(acq_data)
    
    # sort processed acquisition data;
    print('---\n sorting acquisition data...')
    processed_data.sort()
    
    #set the trajectory and compute the dcf
    processed_data = set_grpe_trajectory(processed_data)
    traj = np.transpose(get_grpe_trajectory(processed_data))

    print('---\n computing density weights...')
    dcf = calc_dcf(traj)
    processed_data = set_densitycompensation_as_userfloat(processed_data, dcf)

    # compute coil sensitivity maps
    print('---\n computing coil sensitivity maps...')
    csms = CoilSensitivityData()
    csms.smoothness = 10
    csms.calculate(processed_data)
    
    # create acquisition model based on the acquisition parameters
    print('---\n Setting up Acquisition Model...')

    acq_model = AcquisitionModel()
    acq_model.set_up(processed_data, csms.copy()) #use csm as template image
    acq_model.set_coil_sensitivity_maps(csms)
    
    print('---\n Backward projection ...')
    recon_img = acq_model.backward(processed_data)
    
    if show_plot:
        recon_img.show(title = 'Reconstructed images (magnitude)')

 

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
