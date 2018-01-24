'''Listmode-to-sinograms conversion demo.

Usage:
  listmode_to_sinograms [--help | options] <h_file> <s_file> <t_file>

Arguments:
  h_file  listmode header data file (input)
  s_file  sinogram data file (output)
  t_file  sinogram template data file (input)

Options:
  -p <path>, --path=<path>     path to data files, defaults to data/examples/PET
                               subfolder of SIRF root folder
  -i <int>, --interval=<int>   scanning time interval to convert as string '(a,b)'
                               [default: (0,10)]
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

from ast import literal_eval

from pUtilities import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_path = args['--path']
if data_path is None:
    data_path = petmr_data_path('pet')
prefix = data_path + '/'
h_file = args['<h_file>']
s_file = args['<s_file>']
t_file = args['<t_file>']
interval = literal_eval(args['--interval'])
storage = args['--storage']

def main():

    # select acquisition data storage scheme
    AcquisitionData.set_storage_scheme(storage)

    # create listmode-to-sinograms converter object
    lm2sino = ListmodeToSinograms()

    # set input, output and template files
    lm2sino.set_input(prefix + h_file)
    lm2sino.set_output(s_file)
    lm2sino.set_template(prefix + t_file)

    # set interval
    lm2sino.set_interval(interval[0], interval[1])

    # set flags
    lm2sino.flag_on('store_prompts')
    lm2sino.flag_off('interactive')
    try:
        lm2sino.flag_on('make coffee')
    except error as err:
        print('%s' % err.value)

    # set up the converter
    lm2sino.set_up()

    # convert
    lm2sino.process()

    # get access to the sinograms
    acq_data = lm2sino.get_output()
    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print('acquisition data dimensions: %dx%dx%d' % acq_dim)
    z = acq_dim[0]//2
    show_2D_array('Acquisition data', acq_array[z,:,:])

    # compute randoms
    print('computing randoms, please wait...')
    randoms = lm2sino.compute_randoms()
    rnd_array = randoms.as_array()
    show_2D_array('Randoms', rnd_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
