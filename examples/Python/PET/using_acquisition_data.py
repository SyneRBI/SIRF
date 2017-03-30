'''Acquisition data handling demo.

Usage:
  using_acquisition_data [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -e <engn>, --engine=<engn>  reconstruction engine [default: Stir]
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

from pUtil import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')

def main():

    # output goes to files
    printer = Printer('info.txt', 'warn.txt', 'errr.txt')

    # PET acquisition data to be read from this file
    raw_data_file = existing_filepath(data_path, data_file)
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    z = acq_dim[0]//2

    show_2D_array('Acquisition data', acq_array[z,:,:])

    # clone the acquisition data
    new_acq_data = acq_data.clone()
    # display the cloned data
    acq_array = new_acq_data.as_array()

    show_2D_array('Cloned acquisition data', acq_array[z,:,:])

    # fill the cloned data with the acquisition data multiplied by 10
    new_acq_data.fill(10*acq_array)
    acq_array = new_acq_data.as_array()

    show_2D_array('Scaled acquisition data', acq_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('exception occured: %s' % err.value)
