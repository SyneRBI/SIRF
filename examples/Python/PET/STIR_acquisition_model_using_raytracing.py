'''Demo for setting advanced parameters in sirf.STIR.AcquisitionModelUsingRayTracingMatrix

Usage:
  STIR_acquisition_model_using_raytracing.py

This just sets some parameters and prints settings. Check other demos on how to use it.
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2022 University College London.
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

import sirf.STIR as PET
# create a matrix
m=PET.RayTracingMatrix()
# set some parameters
m.set_do_symmetry_90degrees_min_phi(False)
m.enable_cache(False)
# print information
print(m.get_info())

# same for an acquisition model

# we need to get the matrix itself for the advanced options
am=PET.AcquisitionModelUsingRayTracingMatrix()
am.set_num_tangential_LORs(5)
am.get_matrix().set_restrict_to_cylindrical_FOV(False)
am.get_matrix().set_do_symmetry_90degrees_min_phi(False)
print(am.get_matrix().get_info())
am.get_matrix().set_do_symmetry_180degrees_min_phi(False)
print(am.get_matrix().get_info())
am.get_matrix().set_do_symmetry_swap_s(False)
am.get_matrix().set_do_symmetry_swap_segment(False)
am.get_matrix().set_do_symmetry_shift_z(False)
print(am.get_matrix().get_info())
