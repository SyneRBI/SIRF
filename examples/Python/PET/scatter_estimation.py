'''
Scatter estimation demo.

NOTE: Must be used after running acquisition_sensitivity_from_attenuation.py
(to create attenuation correction factors file) and randoms_from_listmode.py
(to create raw data sinograms file and randoms sinograms file).

Usage:
  scatter_estimation [--help | options]

Options: (defaults are set to work for mMR data processed in the current directory)
  -f <file>, --file=<file>    raw data file [default: sinograms_f1g1d0b0.hs]
  -r <file>, --randoms=<file>  filename with randoms [default: randoms.hs]
  -p <path>, --path=<path>    path to normalization and attenuation files,
                              defaults to data/examples/PET/mMR
  -n <norm>, --norm=<norm>    normalization file [default: norm.n.hdr]
  -a <file>, --attenuation_image=<file>
                              attenuation image file [default: mu_map.hv]
  -A <file>, --attenuation_correction_factors=<file>
                              attenuation correction factors file [default: acf.hs]
  -o <file>, --output=<file>  output prefix for scatter estimates [default: scatter_estimate]
                              ("_#.hs" will be appended, with # the iteration number).
                              Set this to an empty string to prevent output on disk.
  --non-interactive           do not show plots
'''

## CCP SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2019 University of Hull
## Copyright 2020-2021 University College London
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

__version__ = '1.0.1'
from docopt import docopt

args = docopt(__doc__, version=__version__)

# import engine module
import sirf.STIR as PET

from sirf.Utilities import show_2D_array
import PET_plot_functions
#import os


# process command-line options
raw_data_file = args['--file']
randoms_data_file = args['--randoms']
acf_file = args['--attenuation_correction_factors']
data_path = args['--path']
if data_path is None:
    data_path = PET.examples_data_path('PET') + '/mMR'
norm_file = PET.existing_filepath(data_path, args['--norm'])
mu_map_file = PET.existing_filepath(data_path, args['--attenuation_image'])
output_prefix = args['--output']
interactive = not args['--non-interactive']


def main():

    # direct all engine's messages to files
    _ = PET.MessageRedirector('info.txt', 'warn.txt', 'errr.txt')

    PET.AcquisitionData.set_storage_scheme('memory')

    # Create the Scatter Estimator
    # We can use a STIR parameter file like this
    # par_file_path = os.path.join(os.path.dirname(__file__), '..', '..', 'parameter_files')
    # se = PET.ScatterEstimator(PET.existing_filepath(par_file_path, 'scatter_estimation.par'))
    # However, we will just use all defaults here, and set variables below.
    se = PET.ScatterEstimator()

    prompts = PET.AcquisitionData(raw_data_file)
    se.set_input(prompts)
    se.set_attenuation_image(PET.ImageData(mu_map_file))
    if randoms_data_file is None:
        randoms = None
    else:
        randoms = PET.AcquisitionData(randoms_data_file)
        se.set_randoms(randoms)
    if not(norm_file is None):
        se.set_asm(PET.AcquisitionSensitivityModel(norm_file))
    if not(acf_file is None):
        se.set_attenuation_correction_factors(PET.AcquisitionData(acf_file))
    # Could set number of iterations if you want to (we recommend at least 3)
    se.set_num_iterations(1)
    print("number of scatter iterations that will be used: %d" % se.get_num_iterations())
    # Could set number of subsets used by the OSEM reconstruction inside the scatter estimation loop.
    # Here we will set it to 7 (which is in fact the default), which is appropriate for the mMR
    se.set_OSEM_num_subsets(7)
    se.set_output_prefix(output_prefix)
    se.set_up()
    se.process()
    scatter_estimate = se.get_output()

    if not interactive:
        return

    ## show estimated scatter data
    scatter_estimate_as_array = scatter_estimate.as_array()
    show_2D_array('Scatter estimate', scatter_estimate_as_array[0, 0, :, :])

    ## let's draw some profiles to check
    # we will average over all sinograms to reduce noise
    PET_plot_functions.plot_sinogram_profile(prompts, randoms=randoms, scatter=scatter_estimate)

try:
    main()
    print('done')
except PET.error as err:
    print('%s' % err.value)
