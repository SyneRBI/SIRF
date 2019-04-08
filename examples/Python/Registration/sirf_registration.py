'''Regsitration of SIRF images.

Usage:
  registration [--help | options]

Options:
  --eng_ref <eng>              engine for reference image [default: Reg]
  --eng_flo <eng>              engine for floating image [default: Reg]
  --ref <file>                 reference image (default: test.nii.gz)
  --flo <file>                 floating image (default: test2.nii.gz)
  --par <file>                 parameter file (default: niftyreg_aladin.par)
  --algo <algo>                registration algorithm [default: NiftyAladinSym]
  --rmask <file>               mask of reference image
  --fmask <file>               mask of floating image
  --warped <file>              warped image filename [default: output]
  --TM_forward <file>          forward transformation matrix (if rigid/affine)
  --TM_inverse <file>          inverse transformation matrix (if rigid/affine)
  --disp_fwd_4D <file>         4D forward displacement field image
  --def_fwd_4D <file>          4D forward deformation field image
  --disp_inv_4D <file>         4D inverse displacement field image
  --def_inv_4D <file>          4D inverse deformation field image
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 - 2019 University College London.
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

import os

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module
import pReg
exec('import p' + args['--eng_ref'] + ' as eng_ref')
exec('import p' + args['--eng_flo'] + ' as eng_flo')

# process command-line options
ref_file = args['--ref']
flo_file = args['--flo']
par_file = args['--par']
algo = args['--algo']
rmask_file = args['--rmask']
fmask_file = args['--fmask']

# if using the default for any, need to get the examples folder
if (ref_file is None or flo_file is None or par_file is None): 
  SIRF_PATH = os.environ.get('SIRF_PATH')
  if SIRF_PATH is not None:
    examples_path = SIRF_PATH + '/data/examples/Registration'
  else:
    errorMsg = 'You need to set the SIRF_PATH environment variable to allow finding the raw data.'
    raise error(errorMsg)

# reference
if ref_file is None:
    ref_file = examples_path + "/test.nii.gz"

# floating
if flo_file is None:
    flo_file = examples_path + "/test2.nii.gz"

# parameter file
if par_file is None:
    par_file = examples_path + "/paramFiles/niftyreg_aladin.par"


def main():

    # Open reference and floating images
    print('Engine for reference image: ' + args['--eng_ref'])
    print('Engine for floating image: ' + args['--eng_flo'])
    print('Reference image: ' + ref_file)
    print('Floating image: ' + flo_file)
    print('Parameter file: ' + par_file)
    print('Registration algorithm: ' + algo)

    ref = eng_ref.ImageData(ref_file)
    flo = eng_flo.ImageData(flo_file)

    # Dynamically create registration algorithm
    algorithm = getattr(pReg, algo)
    reg = algorithm()
    reg.set_reference_image(ref)
    reg.set_floating_image(flo)
    reg.set_parameter_file(par_file)

    # If masks have been requested, enter them
    if rmask_file:
      rmask = eng_ref.ImageData(rmask_file)
      reg.set_reference_mask(rmask)
    if fmask_file:
      fmask = eng_flo.ImageData(fmask_file)
      reg.set_floating_mask(fmask)

    reg.process()

    # Output
    print("saving to file: " + args['--warped'])
    reg.get_output().write(args['--warped'])

    # TMs
    if args['--TM_forward'] is not False:
        reg.get_transformation_matrix_forward().write(args['--TM_forward'])
    if args['--TM_inverse'] is not False:
        reg.get_transformation_matrix_inverse().write(args['--TM_inverse'])

    # Disp fields
    if args['--disp_fwd_4D'] is not False:
        reg.get_displacement_field_forward().write(args['--disp_fwd_4D'])
    if args['--disp_inv_4D'] is not False:
        reg.get_displacement_field_inverse().write(args['--disp_inv_4D'])

    # Def fields
    if args['--def_fwd_4D'] is not False:
        reg.get_deformation_field_forward().write(args['--def_fwd_4D'])
    if args['--def_inv_4D'] is not False:
        reg.get_deformation_field_inverse().write(args['--def_inv_4D'])


try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
