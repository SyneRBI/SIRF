'''Acquisition data from scanner info demo.

Usage:
  acquisition_data_from_scanner_info [--help | options]

Options:
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
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

from pUtilities import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

storage = args['--storage']

def main():

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # create acquisition data from scanner parameters
    acq_data = AcquisitionData('Siemens_mMR')
    # set all values to 1.0
    acq_data.fill(1.0)

    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print('acquisition data dimensions: %dx%dx%d' % acq_dim)

    z = acq_dim[0]//2
    show_2D_array('Acquisition data', acq_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
