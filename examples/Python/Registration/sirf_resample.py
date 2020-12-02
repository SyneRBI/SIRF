'''Resampling of SIRF images.

Usage:
  resample [--help | options]

Options:
  --eng_ref <eng>              engine for reference image [default: Reg]
  --eng_flo <eng>              engine for floating image [default: Reg]
  --ref <file>                 reference image (default: test.nii.gz)
  --flo <file>                 floating image (default: test2.nii.gz)
  --algo <algo>                resampling algorithm [default: NiftyResample]
  --output <file>              output image filename [default: output]
  --intrp <intrp>              interpolation order, defaults to cubic [default: 3]
  --trans_filenames ...        transformation filenames, (with quotations): "filename1,filename2,filename3"
  --trans_types ...            transformation types, e.g. (with quotations): "AffineTransformation,NiftiImageData3DDeformation,NiftiImageData3DDisplacement"
  --pad <pad>                  Padding value
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 - 2019 University College London.
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

import os

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# import engine module
import sirf.Reg
exec('import p' + args['--eng_ref'] + ' as eng_ref')
exec('import p' + args['--eng_flo'] + ' as eng_flo')

# process command-line options
ref_file = args['--ref']
flo_file = args['--flo']
algo = args['--algo']
pad = args['--pad']

# if using the default for any, need to get the examples folder
if (ref_file or flo_file) is None: 
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

# parse the transformations
trans_filenames_str = args['--trans_filenames']
trans_types_str = args['--trans_types']
trans_filenames = list()
trans_types = list()
if trans_filenames_str:
  trans_filenames = trans_filenames_str.split(",")
if trans_types_str:
  trans_types = trans_types_str.split(",")
if len(trans_filenames) != len(trans_types):
  raise ValueError('Expect as many transformation filenames as types.')

def main():

    # Open reference and floating images
    print('\nEngine for reference image: ' + args['--eng_ref'])
    print('Engine for floating image: ' + args['--eng_flo'])
    print('Reference image: ' + ref_file)
    print('Floating image: ' + flo_file + '\n')

    ref = eng_ref.ImageData(ref_file)
    flo = eng_flo.ImageData(flo_file)

    # Dynamically create resample algorithm. With inline code, you can do e.g. res = sirf.Reg.NiftyResample()
    algorithm = getattr(sirf.Reg, algo)
    res = algorithm()
    # Set the image we want to resample
    res.set_reference_image(ref)
    # the floating image is set so we know the domain of the resampled image. This can be ref.
    res.set_floating_image(flo)
    # 0 is nearest neighbour, 1 is linear, 3 is cubic, 4 is sinc
    res.set_interpolation_type(int(args['--intrp']))

    # create and add each transformation
    for i in range(len(trans_filenames)):
      print('Transformation ' + str(i) + ' filename: ' + trans_filenames[i])
      print('Transformation ' + str(i) + ' type: ' + trans_types[i] + '\n')
      trans_class = getattr(sirf.Reg, trans_types[i])
      trans = trans_class(trans_filenames[i])
      res.add_transformation(trans)

    # If padding value has been set
    if pad is not None:
      res.set_padding_value(pad)

    # Resample
    res.process()

    # Output
    print("saving to file: " + args['--output'])
    res.get_output().write(args['--output'])


try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
