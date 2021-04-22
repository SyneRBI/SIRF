"""Get multiplicative sinograms from normalisation and/or attenuation.

Usage:
  osem_reconstruction [--help | options]

Options:
  -p <path>, --path=<path>      path to data files, defaults to data/examples/PET/mMR
                                subfolder of SIRF root folder
  -S <file>, --sino=<file>      template sinogram [default: mMR_template_span11_small.hs]
  -A <attn>, --attn=<attn>      attenuation image file [default: mu_map.hv]
  -N <norm>, --norm=<norm>      ECAT8 bin normalisation file [default: norm.n.hdr]
  -O <outp>, --outp=<outp>      output file [default: multiplicative]
  -T <file>, --trans=<file>     transform for attn image
  -t <str>, --trans_type=<str>  transform type (tm, disp, def) [default: tm]
"""

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2020 University College London.
#
# This is software developed for the Collaborative Computational
# Project in Synergistic Reconstruction for Biomedical Imaging
# (formerly CCP PETMR)
# (http://www.ccpsynerbi.ac.uk/).
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

from docopt import docopt
from os import path
import sirf.STIR as pet
import sirf.Reg as reg
from sirf.Utilities import error, show_3D_array, examples_data_path, existing_filepath

__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)


def check_file_exists(filename):
    """Check file exists, else throw error."""
    if not path.isfile(filename):
        raise error('File not found: %s' % filename)


# process command-line options
data_path = args['--path']
if data_path is None:
    # default to data/examples/PET/mMR
    # Note: seem to need / even on Windows
    #data_path = os.path.join(examples_data_path('PET'), 'mMR')
    data_path = examples_data_path('PET') + '/mMR'
print('Finding files in %s' % data_path)

# Sinogram. if sino not found, get the one in the example data
sino_file = existing_filepath(data_path, args['--sino'])

# Attenuation - image
attn_im_file = existing_filepath(data_path, args['--attn'])

# Norm - ECAT8
norm_e8_file = existing_filepath(data_path, args['--norm'])

# Attn transformation
trans = args['--trans']
if trans:
    check_file_exists(trans)
trans_type = args['--trans_type']

# Output file
outp_file = args['--outp']


def resample_attn_image(image):
    """Resample the attenuation image."""
    if trans_type == 'tm':
        transformation = reg.AffineTransformation(trans)
    elif trans_type == 'disp':
        transformation = reg.NiftiImageData3DDisplacement(trans)
    elif trans_type == 'def':
        transformation = reg.NiftiImageData3DDeformation(trans)
    else:
        raise ValueError("Unknown transformation type.")

    resampler = reg.NiftyResampler()
    resampler.set_reference_image(image)
    resampler.set_floating_image(image)
    resampler.set_interpolation_type_to_linear()
    resampler.set_padding_value(0.0)
    resampler.add_transformation(transformation)
    return resampler.forward(image)


def main():
    """Do main."""
    # Acq model and template sino
    acq_model = pet.AcquisitionModelUsingRayTracingMatrix()
    acq_data = pet.AcquisitionData(sino_file)

    # If norm is present
    asm_norm = None
    if norm_e8_file:
        # create acquisition sensitivity model from ECAT8 normalisation data
        asm_norm = pet.AcquisitionSensitivityModel(norm_e8_file)

    # If attenuation is present
    asm_attn = None
    if attn_im_file:
        attn_image = pet.ImageData(attn_im_file)
        if trans:
            attn_image = resample_attn_image(attn_image)
        asm_attn = pet.AcquisitionSensitivityModel(attn_image, acq_model)
        # temporary fix pending attenuation offset fix in STIR:
        # converting attenuation into 'bin efficiency'
        asm_attn.set_up(acq_data)
        bin_eff = pet.AcquisitionData(acq_data)
        bin_eff.fill(1.0)
        print('applying attenuation (please wait, may take a while)...')
        asm_attn.unnormalise(bin_eff)
        asm_attn = pet.AcquisitionSensitivityModel(bin_eff)

    # Get ASM dependent on attn and/or norm
    if asm_norm and asm_attn:
        print("AcquisitionSensitivityModel contains norm and attenuation...")
        asm = pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
    elif asm_norm:
        print("AcquisitionSensitivityModel contains norm...")
        asm = asm_norm
    elif asm_attn:
        print("AcquisitionSensitivityModel contains attenuation...")
        asm = asm_attn
    else:
        raise ValueError("Need norm and/or attn")

    # only need to project again if normalisation is added
    # (since attenuation has already been projected)
    if asm_norm:
        asm_attn.set_up(acq_data)
        bin_eff = pet.AcquisitionData(acq_data)
        bin_eff.fill(1.0)
        print('getting sinograms for multiplicative factors...')
        asm.set_up(acq_data)
        asm.unnormalise(bin_eff)

    print('writing multiplicative sinogram: ' + outp_file)
    bin_eff.write(outp_file)


if __name__ == "__main__":
    main()
