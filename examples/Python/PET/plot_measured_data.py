'''
Plot prompts/randoms/scatter demo

Usage:
  plot_measured_data [--help | options]

Options:
  -f <file>, --file=<file>    raw data file [default: sinograms_f1g1d0b0.hs]
  -r <file>, --randoms=<file> filename with randoms [default: None]
  -s <file>, --scatter=<file> scatter data file [default: None]
  --non-interactive           do not show plots
'''

## CCP SyneRBI Synergistic Image Reconstruction Framework (SIRF)
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

__version__ = '1.0.0'
from docopt import docopt

args = docopt(__doc__, version=__version__)
if args['--non-interactive']:
    exit()

prompts = args["--file"]
scatter = args["--scatter"]
randoms = args["--randoms"]

import PET_plot_functions

PET_plot_functions.plot_sinogram_profile(prompts, randoms=randoms, scatter=scatter)
