# CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2018 - 2019 University College London
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

# Imports
import os
import sys
import time
import numpy

import sirf.STIR as pet
import sirf.Gadgetron as mr
import sirf.Reg as reg
from sirf.Utilities import error

# Paths
SIRF_PATH = os.environ.get('SIRF_PATH')


# NiftiImageData
def try_stirtonifti():
    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Starting STIR to Nifti test...\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)

    # Input filenames
    nifti_filename = SIRF_PATH + '/data/examples/Registration/test2.nii.gz';

    # Load the image as a NiftiImageData3D
    image_nifti = reg.NiftiImageData3D(nifti_filename);

    # Read as STIRImageData, convert to NiftiImageData3D and save to file
    image_stir = pet.ImageData(nifti_filename);
    image_nifti_from_stir = reg.NiftiImageData3D(image_stir);
    image_nifti_from_stir.write('results/stir_to_nifti.nii',image_nifti.get_original_datatype());

    # Compare the two
    if image_nifti != image_nifti_from_stir:
        raise AssertionError("Conversion from STIR to Nifti failed.")

    # Resample and then check that voxel values match
    resample = reg.NiftyResample()
    resample.set_floating_image(image_stir) 
    resample.set_reference_image(image_nifti) 
    resample.set_interpolation_type_to_nearest_neighbour()
    resample.process()

    # as_array() of both original images should match
    if not numpy.array_equal(image_nifti.as_array(),resample.get_output().as_array()):
        raise AssertionError("as_array() of sirf.Reg.NiftiImageData and resampled sirf.STIR.ImageData are different.")

    time.sleep(0.5)
    sys.stderr.write('\n# --------------------------------------------------------------------------------- #\n')
    sys.stderr.write('#                             Finished STIR to Nifti test.\n')
    sys.stderr.write('# --------------------------------------------------------------------------------- #\n')
    time.sleep(0.5)


def test():
    try_stirtonifti()


if __name__ == "__main__":
    try:
        test()
    except:
        raise error("Error encountered.")
