'''gpuCgSense demo

Tested on https://sourceforge.net/projects/gadgetron/files/testdata/ismrmrd/golden_angle.h5

Usage:
  golden_radial_mode2_gpusense_cg.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
## Copyright 2015 - 2017 University College London.
## Copyright 2015 - 2017 Physikalisch-Technische Bundesanstalt.
##
## This is software developed for the Collaborative Computational
## Project in Positron Emission Tomography and Magnetic Resonance imaging
## (http://www.ccppetmr.ac.uk/).
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

import time

# import SIRF utilities
from sirf.Utilities import examples_data_path, existing_filepath, error
# import MR engine types
from sirf.Gadgetron import AcquisitionData, Reconstructor

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')

def main():

    # locate the input data
    input_file = existing_filepath(data_path, data_file)
    acq_data = AcquisitionData(input_file)
    
    # create reconstruction object
    recon = Reconstructor([ \
        'NoiseAdjustGadget', \
        'PCACoilGadget', \
        'CoilReductionGadget(coils_out=16)', \
        'gpuRadialSensePrepGadget(' + \
            'mode=2,' + \
            'profiles_per_frame=16,' + \
            'rotations_per_reconstruction=16,' + \
            'buffer_frames_per_rotation=16,' + \
            'buffer_length_in_rotations=2,' + \
            'reconstruction_os_factor_x=1.5,' + \
            'reconstruction_os_factor_y=1.5' + \
        ')', \
        'slice0:gpuCgSenseGadget(' + \
            'number_of_iterations=10,' + \
            'oversampling_factor=1.25,' + \
            'output_convergence=true' + \
        ')', \
        'slice1:gpuCgSenseGadget(' + \
            'sliceno=1,' + \
            'number_of_iterations=10,' + \
            'oversampling_factor=1.25,' + \
            'output_convergence=true' + \
        ')', \
        'slice2:gpuCgSenseGadget(' + \
            'sliceno=2,' + \
            'number_of_iterations=10,' + \
            'oversampling_factor=1.25,' + \
            'output_convergence=true' + \
        ')', \
        'ExtractGadget', 'AutoScaleGadget'])

    # provide raw k-space data as input
    recon.set_input(acq_data)

    # perform reconstruction
    recon.process()

    # retrieve reconstructed image data
    image_data = recon.get_output()

    image_data.show(title = 'Images')

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
