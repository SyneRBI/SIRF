'''Acquisition data from scanner info demo.

Usage:
  acquisition_data_from_scanner_info [--help | options]

Options:
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 Rutherford Appleton Laboratory STFC
## Copyright 2018 University College London.
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
#exec('from p' + args['--engine'] + ' import *')
exec('from sirf.' + args['--engine'] + ' import *')

storage = args['--storage']

def main():

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # create acquisition data from scanner parameters
    acq_data = AcquisitionData('Siemens_mMR')
    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print('acquisition data dimensions (maximum resolution): %dx%dx%dx%d' % acq_dim)

    # create acquisition data from scanner parameters but with axial compression etc
    acq_data = AcquisitionData('Siemens_mMR',span=11, view_mash_factor=2)
    # set all values to 1.0
    acq_data.fill(1.0)
    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print('acquisition data dimensions (span 11, view mashing 2): %dx%dx%dx%d' % acq_dim)

    # By default, the Siemens mMR uses acquisition settings corresponding to the following constructor
    acq_data = AcquisitionData('Siemens_mMR',span=11, max_ring_diff=60, view_mash_factor=1)

    # write the acquisition data to a file (commented out for this demo)
    # acq_data.write('example_mMR_ones.hs')
try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
