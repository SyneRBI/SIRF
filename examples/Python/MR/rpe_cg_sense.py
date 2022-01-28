'''
Example of an iterative reconstruciton with radial phase encoding (RPE) data.

Upper-level demo that illustrates the computation of how to use a non-cartesian
radial phase-encoding acquisition model to reconstruct data iteratively and
without the use of any k-space density weighting.

Usage:
  rpe_recon.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: 3D_RPE_Lowres.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
  -n <bool>, --non-cart=<bool> run recon iff non-cartesian code was compiled
                              [default: False]
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2021 Physikalisch-Technische Bundesanstalt (PTB)
## Copyright 2015 - 2021 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2021 University College London.
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

run_recon = str(args['--non-cart']) == 'True'

import numpy

# define symmetrical operator for cg-optimisation
def EhE(E, image ):
    return E.backward( E.forward(image) )

def main():

    # locate the k-space raw data file
    input_file = existing_filepath(data_path, data_file)
    
    # acquisition data will be read from an HDF file input_file
    # AcquisitionData.set_storage_scheme('memory')
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

    if run_recon is True:
    
        print('---\n computing coil sensitivity maps...')
        csms = CoilSensitivityData()
        csms.smoothness = 10
        csms.calculate(processed_data)
        
        # create acquisition model based on the acquisition parameters
        print('---\n Setting up Acquisition Model...')
    
        #set up the acquisition model
        E = AcquisitionModel()
        E.set_up(processed_data, csms.copy())
        E.set_coil_sensitivity_maps(csms)
    
        print('---\n Backward projection ...')
        # this is our first residual
        recon_img = E.backward(processed_data)
        recon_img.fill(0+0j) # for some reason you need to start with this set to zero
 
        # now copy the pseudo-code from wikipedia for cg optimisation
        x = recon_img
        y = processed_data

        # this is our first residual
        r = E.backward( y ) - EhE(E,x)

        # this is our cost function at the start
        rr = r.norm() ** 2
        rr0 = rr

        # initialize p
        p = r
        
        # define optimisation parameters
        num_iter = 10
        sufficiently_small = 1e-7
        
        print('Cost for k = 0: '  + str( rr/ rr0) )
        
        for k in range(num_iter):

            Ap = EhE(E, p )

            alpha = rr / Ap.dot(p)

            x = x + alpha * p

            r = r - alpha * Ap

            beta  = r.norm()**2 / rr
            rr = r.norm()**2

            p = r + beta * p

            relative_residual = numpy.sqrt(rr/rr0)

            print('Cost at step  {} = {}'.format(k+1, relative_residual))
            
            if( relative_residual  < sufficiently_small ):
                print('We achieved our desired accuracy. Stopping iterative reconstruction')
                break

            if k is num_iter-1:
                print('Reached maximum number of iterations. Stopping reconstruction.')

    else:
        print('---\n Skipping non-cartesian code...')

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)

