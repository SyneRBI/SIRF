'''
Lower-level interface demo that illustrates creating and running a chain
of gadgets - shortest version.

Usage:
  fully_sampled_recon_single_chain_short.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  images output file
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2019 University College London.
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
from pGadgetron import *

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
output_file = args['--output']
show_plot = not args['--non-interactive']

def main():

    # locate the input data
    input_file = existing_filepath(data_path, data_file)
    acq_data = AcquisitionData(input_file)

    # create reconstruction object
    recon = Reconstructor(['RemoveROOversamplingGadget', \
        'SimpleReconGadgetSet'])
    # reconstruct images
    image_data = recon.reconstruct(acq_data)
    # show reconstructed images
    if show_plot:
        image_data.show(title='Reconstructed images (magnitude)')

    if output_file is not None:
        # write images to output_file
        print('writing to %s' % output_file)
        image_data.write(output_file)

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
