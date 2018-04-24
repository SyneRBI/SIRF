'''FBP2D reconstruction demo

Usage:
  fbp2d_reconstruction [--help | options]

Options:
  -d <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -e <engn>, --engn=<engn>    reconstruction engine [default: STIR]
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

# import engine module
exec('from p' + args['--engn'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
raw_data_file = existing_filepath(data_path, data_file)

def main():

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = MessageRedirector()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    recon = FBP2DReconstructor()
    recon.set_input(acq_data)
    recon.reconstruct()
    image = recon.get_output()

    #image_array = image.as_array()
    # interactively display the reconstructed image
    image.show()

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
