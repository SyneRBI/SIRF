'''
Lower-level demo, 2-chain GRAPPA reconstruction of undersampled data.

Usage:
  undersampled_recon_chain.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -g, --gfactors              compute Gfactors
  -o <file>, --output=<file>  images output file
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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

output_file = args['--output']
get_gfactors = args['--gfactors']

# import engine module
from pGadgetron import *

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    # acquisitions will be read from this HDF file
    input_file = existing_filepath(data_path, args['--file'])
    input_data = AcquisitionData(input_file)

    # pre-process acquisitions
    print('---\n pre-processing acquisitions...')
    preprocessed_data = input_data.process(['NoiseAdjustGadget', \
         'AsymmetricEchoAdjustROGadget', 'RemoveROOversamplingGadget'])

    # set up reconstruction chain
    recon = Reconstructor([\
         'AcquisitionAccumulateTriggerGadget', 'BucketToBufferGadget' \
         '(N_dimension=contrast,S_dimension=average,split_slices=false)', \
         'GenericReconCartesianReferencePrepGadget', \
         'GRAPPA:GenericReconCartesianGrappaGadget', \
         'GenericReconFieldOfViewAdjustmentGadget', \
         'GenericReconImageArrayScalingGadget', 'ImageArraySplitGadget'])
    # change a property of the gadget labelled by 'GRAPPA'
    recon.set_gadget_property('GRAPPA', 'send_out_gfactor', get_gfactors)
    recon.set_input(preprocessed_data)
    # reconstruct
    print('---\n reconstructing...')
    recon.process()
    output = recon.get_output()

    # show images
    image_array = output.as_array()
    title = 'Reconstructed images (magnitude)'
    show_3D_array(abs(image_array), suptitle = title)

    if output_file is not None:
        # write images to a new group in args.output
        # named after the current date and time
        time_str = time.asctime()
        print('writing to %s' % output_file)
        output.write(output_file, time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
