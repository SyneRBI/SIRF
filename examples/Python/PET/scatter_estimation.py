'''Scatter estimation demo: Executes the ScatterRun example from
STIR.

Usage:
  scatter_simulation [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: Utahscat600k_ca_seg4.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -a <addv>, --addv=<addv>    additive term value [default: 0]
  -b <back>, --back=<back>    background term value [default: 0]
  -n <norm>, --norm=<norm>    normalization value [default: 1]
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]

There is an interactive demo with much more documentation on this process.
You probably want to check that instead.
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2019 University of Hull
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
exec('from sirf.' + args['--engine'] + ' import *')

from pUtilities import show_2D_array

# import engine module
exec('from sirf.' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
addv = float(args['--addv'])
back = float(args['--back'])
beff = 1 / float(args['--norm'])
output_file = args['--output']


def main():

    # TODO: properly set the path of Scatter Estimation parameter file.
    example_path = os.environ.get('STIR_PATH') + "/recon_test_pack/ScatterRun/"
    os.chdir(example_path)

    # Create the Single Scatter Simulation model
    se = ScatterEstimator('ScatterEstimation.par')

    scatter_estimate = se.process()

    # show simulated scatter data
    scatter_estimate_as_array = scatter_estimate.as_array()
    show_2D_array('Scatter estimation', scatter_estimate_as_array[0, 0, :, :])


try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)