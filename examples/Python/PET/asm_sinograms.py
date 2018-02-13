'''Acquisition sensitivity model using bin efficiencies sinograms demo.

Usage:
  acquisition_sensitivity_model [--help | options]

Options:
  -f <file>, --file=<file>     raw data file [default: my_forward_projection.hs]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
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

import math

from pUtilities import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
storage = args['--storage']

def main():

    # output goes to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # PET acquisition data to be read from this file
    raw_data_file = existing_filepath(data_path, data_file)
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # copy the acquisition data into a Python array and display it
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    z = acq_dim[0]//2
    show_2D_array('Acquisition data', acq_array[z,:,:])

    # create bin efficiencies sinograms
    bin_eff = acq_data.clone()
    bin_eff.fill(2.0)
    bin_eff_arr = bin_eff.as_array()
    bin_eff_arr[:,10:50,:] = 0
    show_2D_array('Bin efficiencies', bin_eff_arr[z,:,:])
    bin_eff.fill(bin_eff_arr)

    # create acquisition sensitivity model from bin efficiencies
    asm = AcquisitionSensitivityModel(bin_eff)

    # apply normalization to acquisition data
    ad = acq_data.clone()
    asm.unnormalise(ad)
    ad_array = ad.as_array()
    show_2D_array('Normalized acquisition data', ad_array[z,:,:])

    # create another bin efficiencies sinograms
    bin_eff_arr[:,10:50,:] = 2.0
    bin_eff_arr[:,60:80,:] = 0
    show_2D_array('Another bin efficiencies', bin_eff_arr[z,:,:])
    bin_eff2 = acq_data.clone()
    bin_eff2.fill(bin_eff_arr)

    # create another acquisition sensitivity model from bin efficiencies
    asm2 = AcquisitionSensitivityModel(bin_eff2)

    # chain the two models
    asm12 = AcquisitionSensitivityModel(asm, asm2)
    asm12.set_up(acq_data)

    # apply the chain of models to acquisition data
    ad = acq_data.clone()
    asm12.unnormalise(ad)
    ad_array = ad.as_array()
    show_2D_array('Chain-normalized acquisition data', ad_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
