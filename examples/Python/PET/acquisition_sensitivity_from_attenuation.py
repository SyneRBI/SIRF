'''Acquisition sensitivity model using attenuation image.

Usage:
  acquisition_sensitivity_from_attenuation [--help | options]

Options:
  -t <temp>, --temp=<temp>     raw data template [default: mMR_template_span11_small.hs]
  -a <attn>, --attn=<attn>     attenuation image file file [default: mu_map.hv]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET/mMR
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
temp_file = args['--temp']
attn_file = args['--attn']
data_path = args['--path']
if data_path is None:
    # default to data/examples/PET/mMR
    # Note: seem to need / even on Windows
    #data_path = os.path.join(examples_data_path('PET'), 'mMR')
    data_path = examples_data_path('PET') + '/mMR'
temp_file = existing_filepath(data_path, temp_file)
attn_file = existing_filepath(data_path, attn_file)
storage = args['--storage']

def main():

    # direct all engine's messages to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # obtain an acquisition data template
    template = AcquisitionData(temp_file)

    # create uniform acquisition data from template
    print('creating uniform acquisition data...')
    acq_data = AcquisitionData(template)
    acq_data.fill(1.0)

    # read attenuation image
    attn_image = ImageData(attn_file)
    attn_image_as_array = attn_image.as_array()
    z = attn_image_as_array.shape[0]//2
    show_2D_array('Attenuation image', attn_image_as_array[z,:,:])

    # create acquisition model
    am = AcquisitionModelUsingRayTracingMatrix()
    am.set_up(template, attn_image)

    # create acquisition sensitivity model from attenuation image
    print('creating acquisition sensitivity model...')
    asm = AcquisitionSensitivityModel(attn_image, am)
    asm.set_up(template)
    am.set_acquisition_sensitivity(asm)
##    print('projecting (please wait, may take a while)...')
##    simulated_data = am.forward(attn_image)

    # apply attenuation to the uniform acquisition data to obtain
    # 'bin efficiencies'
    print('applying attenuation (please wait, may take a while)...')
    asm.unnormalise(acq_data)

    # show 'bin efficiencies'
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    z = acq_dim[0]//2
    show_2D_array('Bin efficiencies', acq_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)

