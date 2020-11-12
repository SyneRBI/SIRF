'''
Scatter estimation demo: Executes the ScatterRun example from STIR.

Usage:
  scatter_estimation [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: sinospan11_f1g1d0b0.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -r <file>, --randoms=<file>  filename with randoms [default: MLrandomsspan11_f1.hs]
  -n <norm>, --norm=<norm>    normalization file [default: norm.n.hdr]
  -a <file>, --attenuation_image=<file>
                              attenuation image file [default: mu.hv]
  -A <file>, --attenuation_correction_factors=<file>
                              attenuation correction factors file [default: acf.hs]
  -o <file>, --output=<file>  output prefix for scatter estimates [default: scatter_estimate]
                              ("_#.hs" will be appended, with # the iteration number).
                              Set this to an empty string to prevent output on disk.

'''

## CCP SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2019 University of Hull
## Copyright 2020 University College London
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

# import engine module
import sirf.STIR as PET

from sirf.Utilities import show_2D_array

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = PET.examples_data_path('PET')
raw_data_file = PET.existing_filepath(data_path, data_file)
randoms_data_file = args['--randoms']
if not(randoms_data_file is None):
    randoms_data_file = PET.existing_filepath(data_path, randoms_data_file)
norm_file = args['--norm']
if not(norm_file is None):
    norm_file = PET.existing_filepath(data_path, norm_file)
acf_file = args['--attenuation_correction_factors']
if not(acf_file is None):
    acf_file = PET.existing_filepath(data_path, acf_file)
mu_map_file = args['--attenuation_image']
if not(mu_map_file is None):
    mu_map_file = PET.existing_filepath(data_path, mu_map_file)
output_prefix = args['--output']


def main():
    import os
    PET.AcquisitionData.set_storage_scheme('memory')

    # TODO: properly set the path of Scatter Estimation parameter file.
    par_file_path = os.path.join(os.path.dirname(__file__), '..', '..', 'parameter_files')
    # Create the Scatter Estimator
    se = PET.ScatterEstimator(PET.existing_filepath(par_file_path, 'scatter_estimation.par'))
  # set/change some parameters here
    se.set_input(PET.AcquisitionData(raw_data_file))
    se.set_attenuation_image(PET.ImageData(PET.existing_filepath(data_path, mu_map_file)))
    if not(randoms_data_file is None):
        se.set_randoms(PET.AcquisitionData(randoms_data_file))
    if not(norm_file is None):
        se.set_asm(PET.AcquisitionSensitivityModel(norm_file))
    if not(acf_file is None):
        se.set_attenuation_correction_factors(PET.AcquisitionData(acf_file))
    se.set_num_iterations(2)
    se.set_output_prefix(output_prefix)
    print("number of iterations that will be used: %d" % se.get_num_iterations())
    se.set_up()
    se.process()
    scatter_estimate = se.get_output()
    # show simulated scatter data
    scatter_estimate_as_array = scatter_estimate.as_array()
    show_2D_array('Scatter estimation', scatter_estimate_as_array[0, 0, :, :])


try:
    main()
    print('done')
except PET.error as err:
    print('%s' % err.value)
