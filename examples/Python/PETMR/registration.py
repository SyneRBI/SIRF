"""Registration demo.

Usage:
  registration [--help | options]

Options:
  -R <file>, --ref=<file>        reference image
  -F <file>, --flo=<file>        floating image
  --num_iters=<val>              number of iterations [default: 10]
  -s <nstp>, --steps=<nstp>      number of steepest descent steps [default: 3]
  --optimal                      use locally optimal steepest ascent
  --cpg_downsample=<val>         factor to downsample control point grid
                                 spacing relative to reference image (e.g., 2
                                 will mean cpg spacing will be double that of
                                 reference image) [default: 2]
  --space=<str>                  space in which to perform registration (image
                                 or sinogram) [default: sinogram]
  -o <str>, --out_prefix=<str>   registered file prefix [default: registered]
  -d <str>, --dvf_prefix=<str>   dvf file prefix [default: dvf]

Options when registering in sinogram space:
  --templ_sino=<file>            template sinogram
  --rands=<file>                 randoms sinogram
  --attn=<file>                  attenuation image
  --norm=<file>                  ECAT8 normalisation file

  --num_subsets=<val>            number of subsets for PET projection
                                 [default: 21]
  --gpu                          use GPU (requires NiftyPET projector)
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
cpg_downsample_factor = float(args['--cpg_downsample'])
sino_file = args['--templ_sino']
rands_file = args['--rands']
attn_file = args['--attn']
norm_file = args['--norm']
num_iters = int(args['--num_iters'])
opt = args['--optimal']
steps = int(args['--steps'])
out_prefix = args['--out_prefix']
dvf_prefix = args['--dvf_prefix']
num_subsets = int(args['--num_subsets'])
registration_space = args['--space']
if registration_space not in {'sinogram', 'image'}:
    raise error("Unknown registration space: " + registration_space)
use_gpu = True if args['--gpu'] else False

if opt:
    import scipy.optimize


def check_file_exists(filename):
    """Check file exists. Else, throw error."""
    if not os.path.isfile(filename):
        raise error("File not found: " + filename)


def get_image(filename, required):
    """Get an image from its filename."""
    if not required and filename is None:
        return None
    check_file_exists(filename)
    if registration_space == 'image':
        try:
            im = sirf.Reg.ImageData(filename)
        except:
            im = sirf.STIR.ImageData(filename)
    else:
        im = sirf.STIR.ImageData(filename)
    return im


def get_sinogram(filename, required):
    """Get an sinogram from its filename."""
    if not required and filename is None:
        return None
    check_file_exists(filename)
    return sirf.STIR.AcquisitionData(filename)


def get_cpg_2_dvf_converter(ref):
    """Get CPG 2 DVF converter."""
    im_spacing = ref.get_geometrical_info().get_spacing()
    cpg_spacing = [spacing * cpg_downsample_factor for spacing in im_spacing]
    cpg_2_dvf_converter = sirf.Reg.ControlPointGridToDeformationConverter()
    cpg_2_dvf_converter.set_cpg_spacing(cpg_spacing)
    cpg_2_dvf_converter.set_reference_image(ref)
    return cpg_2_dvf_converter


def get_dvf(ref):
    """Get initial DVF."""
    disp = sirf.Reg.NiftiImageData3DDisplacement()
    ref_nii = sirf.Reg.NiftiImageData3D(ref)
    disp.create_from_3D_image(ref_nii)
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


def grad_alpha_image(
        reference, floating,
        current_estimate,
        cpg_2_dvf_converter,
        resampler, dvf):
    """Update alpha (b-spline CPG) in image space."""
    # Get gradient of alpha
    grad_obj_fn = current_estimate - reference
    # sirf.Reg.ImageData(grad_obj_fn).write(
    #         "grad_obj_fn")

    grad_dvf = resampler.backward(dvf, floating, grad_obj_fn)
    # grad_dvf.write("grad_dvf")
    # sirf.Reg.NiftiImageData3DDisplacement(grad_dvf).write("grad_dvf_as_disp")
    grad_alpha = cpg_2_dvf_converter.backward(grad_dvf)

    # Return grad alpha
    return grad_alpha


def grad_alpha_sino(
        reference, estimate,
        cpg_2_dvf_converter,
        resampler, alpha,
        attn, emission_acq_model,
        attenuation_acq_model, asm_norm, sino):
    """Update alpha (b-spline CPG) in sinogram space."""
    # Convert alpha to dvf
    dvf = cpg_2_dvf_converter.forward(alpha)
    # Get current estimate (need to resample both emission and attn)
    transformed_estimate = resampler.forward(dvf, estimate)

    # Update the emission_acq_model
    # with the updated attenuation image.
    if attn:
        transformed_attn = resampler.forward(dvf, attn)
        update_asm(emission_acq_model, transformed_attn, asm_norm, sino, ref)

    estimated_data = emission_acq_model.direct(transformed_estimate)
    estimated_projection = \
        estimated_data - emission_acq_model.get_constant_term()
    # Beware: below is -ve PoissonLogLikelihood, which is +ve KL.
    # Just bear this in mind when thinking about maximisation or minimsation,
    # especially when doing JRM (image reconstruction and registration at same
    # time).
    grad_obj = -measured_data / estimated_data + 1
    grad_emission_image = emission_acq_model.backward(grad_obj)
    # New backward method
    grad_dvf_em = resampler.backward(dvf, emission, grad_emission_image)
    if attn:
        grad_attenuation_image = -emission_acq_model.backward(
            grad_obj * estimated_projection)
        # New backward but now with attn. image
        grad_dvf_atn = resampler.backward(dvf, attn, grad_attenuation_image)
        grad_dvf = grad_dvf_em + grad_dvf_atn
    else:
        grad_dvf = grad_dvf_em
    grad_alpha = cpg_2_dvf_converter.backward(grad_dvf)

    # Return grad alpha and the current estimate
    return [grad_alpha, transformed_estimate]


# def get_initial_tau(alpha, grad_alpha):
#     """Get initial tau"""
#     lmd_max = 2*grad_alpha_0.norm()/alpha_0.norm()
#     tau = 1/lmd_max
#     return tau


# def update_tau(alpha, alpha_0, grad_alpha, grad_alpha_0, max_step):
#     """Get updated tau."""
#     d_alpha = alpha - alpha_0
#     d_grad = grad_alpha - grad_alpha_0
#     # dg = H di, hence a rough idea about lmd_max is given by
#     lmd_max = 2*d_grad.norm()/d_alpha.norm()
#     # alternative smaller estimate for lmd_max is
#     # lmd_max = -2*dg.dot(di)/di.dot(di)
#     tau = min(tau_0, 1/lmd_max)


def main():
    """Do the main function."""

    # Read input files
    ref = get_image(ref_file, True)
    flo = get_image(flo_file, True)

    if registration_space == 'sinogram':
        template_sino = get_sinogram(sino_file, True)
        attn = get_image(attn_file, False)
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

    if registration_space == 'sinogram':
        # Get acquisition models. We'll need 2 - one for emission (which will
        # contain all corrections: randoms, norms, up-to-date attn), and the
        # other for attenuation (blank just for line integrals).
        emission_acq_model = get_acq_model(ref, template_sino)
        attenuation_acq_model = get_acq_model(ref, template_sino)
        attenuation_acq_model.set_up(template_sino, ref)
        update_asm(acq_model, attn, asm_norm, template_sino, ref)

        sino = acq_model.forward(ref)

    # Initial current estimate
    current_estimate = resampler.forward(dvf, flo)
    sirf.Reg.ImageData(current_estimate).write(out_prefix + "_0")
    dvf.write(dvf_prefix + "_0")

    # Optimisation loop
    for iter in range(num_iters):

        # update alpha (b-spline CPG) in image or sinogram space
        # Image space
        if registration_space == 'image':
            grad_alpha = grad_alpha_image(
                ref, flo, current_estimate,
                cpg_2_dvf_converter, resampler, dvf)
        # Sinogram space
        else:
            grad_alpha = grad_alpha_sino(
                ref, flo, current_estimate, cpg_2_dvf_converter, resampler,
                dvf, attn, emission_acq_model, attenuation_acq_model,
                asm_norm, sino)

        # update alpha
        grad_alpha_max = grad_alpha.as_array().max()
        if grad_alpha_max <= 0:
            raise error("damn.")
        print("Max alpha: " + str(alpha.as_array().max()))
        print("Max grad_alpha: " + str(grad_alpha_max))
        alpha -= grad_alpha / grad_alpha.as_array().max()

        # Get current dvf and estimate
        dvf = cpg_2_dvf_converter.forward(alpha)
        resampler.forward(dvf, flo, current_estimate)

        sirf.Reg.ImageData(current_estimate).write(
            out_prefix + "_" + str(iter+1))
        sirf.Reg.NiftiImageData3DDisplacement(dvf).write(
            dvf_prefix + "_" + str(iter+1))


if __name__ == "__main__":
    main()
