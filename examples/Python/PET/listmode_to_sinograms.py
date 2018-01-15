'''ListmodeToSinogram demo.

Usage:
  lm2sino [--help | options] <h_file> <s_file> <t_file>

Arguments:
  h_file  listmode header data file (input)
  s_file  sinogram data file (output)
  t_file  sinogram template data file (input)

Options:
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

from shutil import copyfile

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
storage = args['--storage']

### quick fix for data path problem in the listmode header
##i = h_file.find('.')
##l_file = h_file[0:i] + '.l'
##copyfile(prefix + l_file, l_file)

def main():

    lm2sino = ListmodeToSinograms()
    lm2sino.set_input(prefix + h_file)
    lm2sino.set_output(s_file)
    lm2sino.set_template(prefix + t_file)
    lm2sino.set_interval(0, 10)
    lm2sino.set_up()

    lm2sino.process()

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
