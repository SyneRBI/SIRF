"""Registration demo.

Usage:
  registration [--help | options]

Options:
  -R <file>, --ref=<file>    reference image
  -F <file>, --flo=<file>    floating image
  --cpg_downsample=<val>     factor to downsample control point grid spacing
                             relative to reference image (e.g., 2 will mean
                             cpg spacing will be double that of reference
                             image) [default: 2]
  --templ_sino=<file>        template sinogram
  --rands=<file>             randoms sinogram
  --attn=<file>              attenuation image
  --norm=<file>              ECAT8 normalisation file
  --num_iters=<val>          number of iterations [default: 10]
  --num_subsets=<val>        number of subsets for PET projection [default: 21]
  --space=<str>              space in which to perform registration (image or
                             sinogram) [default: sinogram]
  --gpu                      use GPU
"""

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2020 University College London.
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

import os
from docopt import docopt
from tqdm import tqdm
from sirf.Utilities import error
import sirf.STIR
import sirf.Reg

__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)


# process command-line options
ref_file = args['--ref']
flo_file = args['--flo']
ref_eng = args['--ref_eng']
flo_eng = args['--flo_eng']
cpg_downsample_factor = float(args['--cpg_downsample'])
sino_file = args['--templ_sino']
rands_file = args['--rands']
attn_file = args['--attn']
norm_file = args['--norm']
num_iters = float(args['--num_iters'])
num_subsets = float(args['--num_subsets'])
registration_space = args['--space']
if registration_space not in {'sinogram', 'image'}:
    raise error("Unknown registration space: " + registration_space)
use_gpu = True if args['--gpu'] else False


def check_file_exists(fname):
    """Check file exists. Else, throw error."""
    if not os.path.isfile(filename):
        raise error("File not found: " + filename)


def get_image(filename):
    """Get an image from its filename."""
    check_file_exists(filename)
    return sirf.STIR.ImageData(filename)


def get_sinogram(filename):
    """Get an sinogram from its filename."""
    check_file_exists(filename)
    return sirf.STIR.AcquisitionData(filename)


def get_cpg_2_dvf_converter(ref):
    """Get CPG 2 DVF converter."""
    cpg_spacing = ref.get_voxel_sizes()[1:4] * cpg_downsample_factor
    cpg_2_dvf_converter = sirf.Reg.ControlPointGridToDeformationConverter()
    cpg_2_dvf_converter.set_cpg_spacing(cpg_spacing)
    cpg_2_dvf_converter.set_reference_image(ref)
    return cpg_2_dvf_converter


def get_dvf(ref):
    """Get initial DVF."""
    disp = sirf.Reg.NiftiImageData3DDisplacement()
    disp.create_from_3D_image(ref)
    disp.fill(0)
    dvf = sirf.Reg.NiftiImageData3DDeformation(disp)
    return dvf


def get_resampler(ref, flo):
    """Get resampler."""
    nr = sirf.Reg.NiftyResample()
    nr.set_reference_image(ref)
    nr.set_floating_image(flo)
    nr.set_interpolation_type_to_linear()
    nr.set_padding_value(0)
    return nr


def get_asm_norm(norm_file):
    if norm_file is None:
        return None
    check_file_exists(norm_file)
    asm_norm = pet.AcquisitionSensitivityModel(norm_file)
    return asm_norm


def get_asm_attn(sino, attn, acq_model):
    """Get attn ASM from sino, attn image and acq model"""
    if attn is None:
        return None
    asm_attn = pet.AcquisitionSensitivityModel(attn, acq_model)
    # temporary fix pending attenuation offset fix in STIR:
    # converting attenuation into 'bin efficiency'
    asm_attn.set_up(sino)
    bin_eff = pet.AcquisitionData(sino)
    bin_eff.fill(1.0)
    asm_attn.unnormalise(bin_eff)
    asm_attn = pet.AcquisitionSensitivityModel(bin_eff)
    return asm_attn


def get_acq_model(ref, sino):
    """Get acquisition model."""
    if not use_gpu:
        acq_model = sirf.Reg.AcquisitionModelUsingRayTracingMatrix()
    else:
        acq_model = sirf.Reg.AcquisitionModelUsingNiftyPET()
        acq_model.set_use_truncation(True)
        acq_model.set_cuda_verbosity(0)

    # Add randoms if desired
    if rands_file:
        rands = get_sinogram(rands_file)
        acq_model.set_background_term(rands)

    return acq_model


def update_asm(acq_model, attn, asm_norm, sino, ref):
    """Update acquisition sensitivity models."""

    # Create attn ASM if necessary
    asm_attn = get_asm_attn(attn)

    asm = None
    if asm_norm and asm_attn:
        asm = pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
    elif asm_norm:
        asm = asm_norm
    elif asm_attn:
        asm = asm_attn
    if asm:
        acq_model.set_acquisition_sensitivity(asm)

    # Set up
    acq_model.set_up(sino, ref)


def update_alpha_image(reference, current_estimate,
                       cpg_2_dvf_converter,
                       resampler, bsplines, alpha):
    """Update alpha (b-spline CPG) in image space."""
    # Convert alpha to dvf
    dvf = cpg_2_dvf_converter.forward(alpha)
    # Get current estimate
    transformed_estimate = resampler.forward(dvf, current_estimate)
    # Get gradient of alpha
    grad_obj_fn = (reference - transformed_estimate) * transformed_estimate
    grad_dvf = resampler.backward(grad_obj_fn)
    grad_alpha = cpg_2_dvf_converter.backward(grad_dvf)
    return grad_alpha


def update_alpha_sino(reference, current_estimate,
                      cpg_2_dvf_converter,
                      resampler, bsplines, alpha,
                      attn, acq_model, asm_norm, sino):
    """Update alpha (b-spline CPG) in sinogram space."""
    # Convert alpha to dvf
    dvf = cpg_2_dvf_converter.forward(alpha)
    # Get current estimate (need to resample both emission and attn)
    transformed_estimate = resampler.forward(dvf, current_estimate)
    transformed_attn = resampler.forward(dvf, attn)

    # Update the acq_model and objective function
    # with the updated attenuation image.
    update_asm(acq_model, transformed_attn, asm_norm, sino, ref)

    estimated_projection = acq_model.forward(transformed_estimate)
    estimated_data = estimated_projection + acq_model.get_background_term()
    grad_obj = -measured_data / estimated_data + 1
    grad_emission_image = acq_model.backward(grad_obj)
    # New backward method
    grad_dvf_em = resampler.backward(grad_emission_image)
    # atn_acq_model - AcqModel for everything incl. attn
    grad_attenuation_image = -acq_model.backward(
        grad_obj * estimated_projection)
    # New backward but now with attn. image
    grad_dvf_atn = resampler.backward(grad_attenuation_image)
    grad_dvf = grad_dvf_em + grad_dvf_atn
    grad_alpha = cpg_2_dvf_converter.backward(grad_dvf)
    return grad_alpha


def main():
    """Do the main function."""

    # Read input files
    ref = get_image(ref_file, ref_eng)
    flo = get_image(flo_file, ref_eng)
    template_sino = get_sinogram(sino_file)
    attn = get_image(attn_file) if attn_file else None
    asm_norm = get_asm_norm(norm_file)

    # Ability to convert between CPG and DVF
    cpg_2_dvf_converter = get_cpg_2_dvf_converter(ref)

    # Create DVF and CPG (called alpha)
    dvf = get_dvf(ref)
    alpha = cpg_2_dvf_converter.backward(dvf)

    # We'll need a niftyreg resampler
    nr = get_resampler(ref, flo)

    # And we'll need the ImageGradientWRTDeformationTimesImage class
    # (which we'll call the resampler)
    resampler = sirf.Reg.ImageGradientWRTDeformationTimesImage()
    resampler.set_resampler(nr)

    # Get likelihood
    acq_model = get_acq_model(ref, template_sino)
    update_asm(acq_model, attn, asm_norm, template_sino, ref)

    # Registered image starts as copy of reference image
    registered_im = ref.copy()

    sino = acq_model.forward(ref)
    simulated_sino = sino.copy()

    # Optimisation loop
    for iter in tqdm(range(num_iters)):

        # update alpha (b-spline CPG) in image or sinogram space
        if regularisation_space == 'image':
            alpha = update_alpha_image(reference, current_estimate,
                                       cpg_2_dvf_converter,
                                       resampler, alpha)
        else:
            alpha = update_alpha_sino(reference, current_estimate,
                                      cpg_2_dvf_converter,
                                      resampler, alpha,
                                      attn, acq_model, asm_norm, sino)

        alpha + step * grad_alpha


if __name__ == "__main__":
    main()
