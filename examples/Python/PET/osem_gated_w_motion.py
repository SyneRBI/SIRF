"""OSEM reconstruction demo for gated data with motion.

We actually use the OSMAPOSL reconstructor in this demo. This reconstructor
implements an Ordered Subsets (OS) version of the One Step Late algorithm (OSL)
from Green et al for Maximum a Posteriori (MAP) maximisation. Here we use it
for Maximum Likelihood (ML) in which case it is equivalent to OSEM.

Usage:
  osem_reconstruction_gated_w_motion [--help | options]

Options:
  -S <file>, --sino=<file>          sinogram pattern, * or % wildcard
                                    (e.g., sino_ms*.hs). Enclose in quotations.
  -a <pattern>, --attn=<pattern>    attenuation pattern, * or % wildcard
                                    (e.g., attn_ms*.hv). Enclose in quotations.
  -R <pattern>, --rand=<pattern>    randoms pattern, * or % wildcard
                                    (e.g., rand_ms*.hs). Enclose in quotations.
  -n <norm>, --norm=<norm>          ECAT8 bin normalization file
  -T <pattern>, --trans=<pattern>   transformation pattern, * or % wildcard
                                    (e.g., tm_ms*.txt). Enclose in quotations.
  -t <str>, --trans_type=<str>      transformation type (tm, disp, def)
                                    [default: tm]
  -s <subs>, --subs=<subs>          number of subsets [default: 21]
  -i <siter>, --subiter=<siter>     number of sub-iterations [default: 2]
  -e <engn>, --engine=<engn>        reconstruction engine [default: STIR]
  -o <outp>, --outp=<outp>          output file prefix [default: recon]
  --initial=<str>                   initial estimate
  --nxny=<nxny>                     image x and y dimension [default: 127]
  --dxdy=<dxdy>                     image x and y spacing
                                    (default: determined by scanner)
  --b_spline_order=<val>            b-spline order for resampler [default: 1]
  --verbosity=<int>                 Verbosity [default: 0]
  --nifti                           save output as nifti
  --gpu                             use gpu
"""

# SyneRBI Synergistic Image Reconstruction Framework (SIRF)
# Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
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
from glob import glob
from docopt import docopt
import sirf.Reg as reg


__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)

# import engine module
if args['--engine'] == 'STIR':
    import sirf.STIR as pet

pet.AcquisitionData.set_storage_scheme('memory')

# Multiple files
sino_pattern = str(args['--sino']).replace('%', '*')
attn_pattern = str(args['--attn']).replace('%', '*')
rand_pattern = str(args['--rand']).replace('%', '*')
trans_pattern = str(args['--trans']).replace('%', '*')
trans_type = str(args['--trans_type'])

if attn_pattern is None:
    attn_pattern = ""
if rand_pattern is None:
    rand_pattern = ""

# Norm
norm_file = None
if args['--norm']:
    norm_file = str(args['--norm'])
    if not os.path.isfile(norm_file):
        raise pet.error("Norm file not found: " + norm_file)

# Initial estimate
initial_estimate = None
if args['--initial']:
    initial_estimate = str(args['--initial'])

# Verbosity
verbosity = int(args['--verbosity'])
pet.set_verbosity(verbosity)
if verbosity == 0:
    msg_red = pet.MessageRedirector(None, None, None)

# Rest of the input options
num_subsets = int(args['--subs'])
num_subiterations = int(args['--subiter'])
nxny = int(args['--nxny'])
outp_file = args['--outp']
b_spline_order = args['--b_spline_order']
nifti = True if args['--nifti'] else False
use_gpu = True if args['--gpu'] else False


def get_resampler(image, ref=None, trans=None):
    """Returns a NiftyResample object for the specified transform and image."""
    if ref is None:
        ref = image
    resampler = reg.NiftyResample()
    resampler.set_reference_image(ref)
    resampler.set_floating_image(image)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    if trans is not None:
        resampler.add_transformation(trans)
    return resampler


def get_asm_attn(sino, attn, acq_model):
    """Get attn ASM from sino, attn image and acq model."""
    asm_attn = pet.AcquisitionSensitivityModel(attn, acq_model)
    # temporary fix pending attenuation offset fix in STIR:
    # converting attenuation into 'bin efficiency'
    asm_attn.set_up(sino)
    bin_eff = pet.AcquisitionData(sino)
    bin_eff.fill(1.0)
    asm_attn.unnormalise(bin_eff)
    asm_attn = pet.AcquisitionSensitivityModel(bin_eff)
    return asm_attn


def parse_input_files():
    """Parse input files."""

    if trans_pattern is None:
        raise AssertionError("--trans missing")
    if sino_pattern is None:
        raise AssertionError("--sino missing")
    trans_files = sorted(glob(trans_pattern))
    sino_files = sorted(glob(sino_pattern))
    attn_files = sorted(glob(attn_pattern))
    rand_files = sorted(glob(rand_pattern))

    num_ms = len(sino_files)
    # Check some sinograms found
    if num_ms == 0:
        raise AssertionError("No sinograms found!")
    # Should have as many trans as sinos
    if num_ms != len(trans_files):
        raise AssertionError("#trans should match #sinos. "
                             "#sinos = " + str(num_ms) +
                             ", #trans = " + str(len(trans_files)))
    # If any rand, check num == num_ms
    if len(rand_files) > 0 and len(rand_files) != num_ms:
        raise AssertionError("#rand should match #sinos. "
                             "#sinos = " + str(num_ms) +
                             ", #rand = " + str(len(rand_files)))

    # For attn, there should be 0, 1 or num_ms images
    if len(attn_files) > 1 and len(attn_files) != num_ms:
        raise AssertionError("#attn should be 0, 1 or #sinos")

    return num_ms, trans_files, sino_files, attn_files, rand_files


def read_trans(num_ms, trans_files):
    """Read transformations."""
    trans = num_ms*[None]
    for i in range(num_ms):
        if trans_type == "disp":
            trans[i] = reg.NiftiImageData3DDisplacement(trans_files[i])
        elif trans_type == "tm":
            trans[i] = reg.NiftiImageData3DDisplacement(
                reg.AffineTransformation(trans_files[i]))
        elif trans_type == "def":
            trans[i] = reg.NiftiImageData3DDisplacement(
                reg.NiftiImageData3DDeformation(trans_files[i]))
        else:
            raise pet.error("Unknown transformation type")
    return trans


def get_initial_estimate(sinos):
    """Get initial estimate."""
    if initial_estimate:
        image = pet.ImageData(initial_estimate)
    else:
        # Create image based on ProjData
        image = sinos[0].create_uniform_image(0.0, (nxny, nxny))
        # If using GPU, need to make sure that image is right size.
        if use_gpu:
            dim = (127, 320, 320)
            spacing = (2.03125, 2.08626, 2.08626)
        # elif non-default spacing desired
        elif args['--dxdy']:
            dim = image.dimensions()
            dxdy = float(args['--dxdy'])
            spacing = (image.voxel_sizes()[0], dxdy, dxdy)
        if use_gpu or args['--dxdy']:
            image.initialise(dim=dim,
                             vsize=spacing)
    return image


def get_attn_images(num_ms, attns, trans, image):
    """Get attn images.

    If none supplied, return None.
    If 1 supplied, resample to each motion state.
    If num_ms supplied return them.
    If GPU projector, resample attn images for correct dimensions.
    """
    resampled_attns = None
    if len(attns) > 0:
        resampled_attns = [0]*num_ms
        # if using GPU, dimensions of attn and recon images have to match
        ref = image if use_gpu else None
        for i in range(num_ms):
            # if we only have 1 attn image, then we need to resample into
            # space of each gate. However, if we have num_ms attn images, then
            # assume they are already in the correct position, so use None as
            # transformation.
            tran = trans[i] if len(attns) == 1 else None
            # If only 1 attn image, then resample that. If we have num_ms attn
            # images, then use each attn image of each frame.
            attn = attns[0] if len(attns) == 1 else attns[i]
            resam = get_resampler(attn, ref=ref, trans=tran)
            resampled_attns[i] = resam.forward(attn)
    return resampled_attns


def get_acq_models(num_ms, resampled_attns, sinos, rands, image):
    """Get acquisition models."""
    print("Setting up acquisition models...")
    if not use_gpu:
        acq_models = num_ms * [pet.AcquisitionModelUsingRayTracingMatrix()]
    else:
        acq_models = num_ms * [pet.AcquisitionModelUsingNiftyPET()]
        for acq_model in acq_models:
            acq_model.set_use_truncation(True)
            acq_model.set_cuda_verbosity(verbosity)

    # If present, create ASM from ECAT8 normalisation data
    asm_norm = None
    if norm_file:
        asm_norm = pet.AcquisitionSensitivityModel(norm_file)

    # Loop over each motion state
    for ind in range(num_ms):
        # Create attn ASM if necessary
        asm_attn = None
        if resampled_attns:
            asm_attn = get_asm_attn(sinos[ind], resampled_attns[ind],
                                    acq_models[ind])

        # Get ASM dependent on attn and/or norm
        asm = None
        if asm_norm and asm_attn:
            if ind == 0:
                print("ASM contains norm and attenuation...")
            asm = pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
        elif asm_norm:
            if ind == 0:
                print("ASM contains norm...")
            asm = asm_norm
        elif asm_attn:
            if ind == 0:
                print("ASM contains attenuation...")
            asm = asm_attn
        if asm:
            acq_models[ind].set_acquisition_sensitivity(asm)

        if len(rands) > 0:
            acq_models[ind].set_background_term(rands[ind])

        # Set up
        acq_models[ind].set_up(sinos[ind], image)
        return acq_models


def set_up_reconstructor(num_ms, sinos, acq_models, trans, image):
    """Set up reconstructor."""
    print("Setting up reconstructor...")

    # define objective function to be maximized as
    # Poisson logarithmic likelihood (with linear model for mean)
    obj_fun = pet.PoissonLogLhLinModMeanGatedProjDataWMotion()

    # Loop over each motion state and set
    for i in range(num_ms):
        obj_fun.add_gate(sinos[i], acq_models[i], trans[i])

    # select Ordered Subsets Maximum A-Posteriori One Step Late as the
    # reconstruction algorithm (since we are not using a penalty, or prior, in
    # this example, we actually run OSEM);
    # this algorithm does not converge to the maximum of the objective function
    # but is used in practice to speed-up calculations
    recon = pet.OSMAPOSLReconstructor()
    recon.set_objective_function(obj_fun)
    recon.set_num_subsets(num_subsets)
    recon.set_num_subiterations(num_subiterations)

    # set up the reconstructor based on a sample image
    # (checks the validity of parameters, sets up objective function
    # and other objects involved in the reconstruction, which involves
    # computing/reading sensitivity image etc etc.)
    print('setting up, please wait...')
    recon.set_up(image)

    # set the initial image estimate
    recon.set_current_estimate(image)
    return recon


def main():

    ###########################################################################
    # Parse input files
    ###########################################################################

    num_ms, trans_files, sino_files, attn_files, rand_files = \
        parse_input_files()

    ###########################################################################
    # Read input
    ###########################################################################

    trans = read_trans(num_ms, trans_files)
    sinos = [pet.AcquisitionData(file) for file in sino_files]
    attns = [pet.ImageData(file) for file in attn_files]
    rands = [pet.AcquisitionData(file) for file in rand_files]

    ###########################################################################
    # Initialise recon image
    ###########################################################################

    image = get_initial_estimate(sinos)

    ###########################################################################
    # Resample attenuation images (if necessary)
    ###########################################################################

    resampled_attns = get_attn_images(num_ms, attns, trans, image)

    ###########################################################################
    # Set up acquisition models
    ###########################################################################

    acq_models = get_acq_models(num_ms, resampled_attns, sinos, rands, image)

    ###########################################################################
    # Set up reconstructor
    ###########################################################################

    recon = set_up_reconstructor(num_ms, sinos, acq_models, trans, image)

    ###########################################################################
    # Reconstruct
    ###########################################################################

    print('reconstructing, please wait...')
    recon.process()
    out = recon.get_output()
    if not args['--nifti']:
        out.write(outp_file)
    else:
        reg.NiftiImageData(out).write(outp_file)


# if anything goes wrong, an exception will be thrown
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except pet.error as err:
    # display error information
    print('%s' % err.value)
