'''Acquisition sensitivity model using ECAT8 bin normalization.

Usage:
  acquisition_sensitivity_from_ecat8 [--help | options]

Options:
  -t <temp>, --temp=<temp>     raw data template [default: mMR_template_span11.hs]
  -n <norm>, --norm=<norm>     ECAT8 bin normalization file [default: norm.n.hdr]
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET/mMR
                               subfolder of SIRF root folder
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
  --non-interactive            do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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

import math

from sirf.Utilities import show_2D_array

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')


# process command-line options
temp_file = args['--temp']
norm_file = args['--norm']
data_path = args['--path']
if data_path is None:
    # default to data/examples/PET/mMR
    # Note: seem to need / even on Windows
    #data_path = os.path.join(examples_data_path('PET'), 'mMR')
    data_path = examples_data_path('PET') + '/mMR'
temp_file = existing_filepath(data_path, temp_file)
norm_file = existing_filepath(data_path, norm_file)
storage = args['--storage']
show_plot = not args['--non-interactive']


def main():

    # direct all engine's messages to files
    msg_red = MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # obtain an acquisition data template
    template = AcquisitionData(temp_file)

    # create a uniform acquisition data from template
    acq_data = AcquisitionData(template)
    acq_data.fill(1.0)

    # create acquisition sensitivity model from ECAT8 normalization data
    asm = AcquisitionSensitivityModel(norm_file)
    asm.set_up(template)

    # apply normalization to the uniform acquisition data to obtain
    # bin efficiencies
    fwd_data = asm.forward(acq_data)

    # show bin efficiencies
    acq_array = fwd_data.as_array()
    acq_dim = acq_array.shape
    z = acq_dim[1]//2
    if show_plot:
        show_2D_array('Bin efficiencies', acq_array[0,z,:,:])


try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    print('%s' % err.value)
