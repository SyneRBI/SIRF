"""Save a sinogram as a Nifti image.

This might be useful for observing with an external viewer.
Geometrical information, such as spacing, is set arbitrarily.


Usage:
  sirf_sino_as_im_nii [--help | options]

Options:
  -s <file>, --sino=<file>  input sinogram
  -i <path>, --img=<path>   output image [default: im]
"""

# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2020 University College London.
#
# This is software developed for the Collaborative Computational
# Project in Positron Emission Tomography and Magnetic Resonance imaging
# (http://www.ccppetmr.ac.uk/).
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#       http://www.apache.org/licenses/LICENSE-2.0
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from sirf.Utilities import error
import sirf.STIR as pet
import sirf.Reg as reg
from docopt import docopt

__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)

# Get filenames
f_sino = args['--sino']
f_im = args['--img']


def main():

    sino = pet.AcquisitionData(f_sino)
    image = sino.convert_to_image()
    reg.ImageData(image).write(f_im)


# if anything goes wrong, an exception will be thrown
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
