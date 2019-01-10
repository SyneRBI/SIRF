'''Demo for synergistic resampling.

Usage:
  synergistic_resample [--help | options]

Options:
  --eng_ref <eng>              engine for reference image [default: Reg]
  --eng_flo <eng>              engine for floating image [default: Reg]
  --eng_reg <eng>              engine for registration/resampling package [default: Reg]
  --ref <file>                 reference image
  --flo <file>                 floating image
  --output <file>              output image filename [default: output]
  --intrp <intrp>              interpolation order, defaults to cubic [default: 3]
  --aff <file>                 affine transformation matrix
  --disp <file>                displacement field image
  --deff <file>                deformation field image
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 Rutherford Appleton Laboratory STFC
## Copyright 2018 University College London.
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
exec('import p' + args['--eng_ref'] + ' as eng_ref')
exec('import p' + args['--eng_flo'] + ' as eng_flo')
exec('import p' + args['--eng_reg'] + ' as eng_reg')

# process command-line options
SIRF_PATH = os.environ.get('SIRF_PATH')
if SIRF_PATH is not None:
  examples_path = SIRF_PATH + '/data/examples/Registration'
else:
  errorMsg = 'You need to set the SIRF_PATH environment variable to allow finding the raw data.'
  raise error(errorMsg)

# reference
ref_file = args['--ref']
if ref_file is None:
  ref_file = examples_path + "/test.nii.gz"

# floating
flo_file = args['--flo']
if flo_file is None:
  flo_file = examples_path + "/test2.nii.gz"

# If no transformations given, use identity
if args['--aff'] is not False:
  trans = eng_reg.AffineTransformation(args['--aff'])
elif args['--disp'] is not False:
  eng_reg.NiftiImageData3DDisplacement(args['--disp'])
elif args['--def'] is not False:
  eng_reg.NiftiImageData3DDisplacement(args['--def'])
else:
  trans = eng_reg.AffineTransformation.get_identity()


def main():

    # Open reference and floating images
    print('Engine for reference image: ' + args['--eng_ref'])
    print('Engine for floating image: ' + args['--eng_flo'])
    print('Engine for resampling: ' + args['--eng_reg'])
    print('Reference image: ' + ref_file)
    print('Floating image: ' + flo_file)

    ref = eng_ref.ImageData(ref_file)
    flo = eng_flo.ImageData(flo_file)

    # Resample
    nr = eng_reg.NiftyResample()
    nr.set_reference_image(ref)
    nr.set_floating_image(flo)
    nr.set_interpolation_type(int(args['--intrp']))
    nr.add_transformation(trans)
    nr.process()

    # Output
    print("saving to file: " + args['--output'])
    nr.get_output().write(args['--output'])


try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
