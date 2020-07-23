"""MCIR for PET.

Usage:
  PET_MCIR [--help | options]

Options:
  -T <pattern>, --trans=<pattern>   transformation pattern, * or % wildcard
                                    (e.g., tm_ms*.txt). Enclose in quotations.
  -t <str>, --trans_type=<str>      transformation type (tm, disp, def)
                                    [default: tm]
  -S <pattern>, --sino=<pattern>    sinogram pattern, * or % wildcard
                                    (e.g., sino_ms*.hs). Enclose in quotations.
  -a <pattern>, --attn=<pattern>    attenuation pattern, * or % wildcard
                                    (e.g., attn_ms*.hv). Enclose in quotations.
  -R <pattern>, --rand=<pattern>    randoms pattern, * or % wildcard
                                    (e.g., rand_ms*.hs). Enclose in quotations.
  -n <norm>, --norm=<norm>          ECAT8 bin normalization file
  -i <int>, --iter=<int>            num iterations [default: 10]
  -r <string>, --reg=<string>       regularisation ("none","FGP_TV", ...)
                                    [default: none]
  -o <outp>, --outp=<outp>          output file prefix [default: recon]
  --nxny=<nxny>                     image x and y dimension [default: 127]
  --dxdy=<dxdy>                     image x and y spacing
                                    (default: determined by scanner)
  -I <str>, --initial=<str>         Initial estimate
  --visualisations                  show visualisations
  --nifti                           save output as nifti
  --gpu                             use GPU projector
  -v <int>, --verbosity=<int>       STIR verbosity [default: 0]
  -s <int>, --save_interval=<int>   save every x iterations [default: 10]
  --descriptive_fname               option to have descriptive filenames
  --update_obj_fn_interval=<int>    frequency to update objective function
                                    [default: 1]
  --sigma=<val>                     set PDHG sigma (default: 1/normK)
  --tau=<val>                       set PDHG tau (default: 1/normK)
                                    Calculating this takes time.
  --alpha=<val>                     regularisation strength (if used)
                                    [default: 0.5]
  --reg_iters=<val>                 Number of iterations for the regularisation
                                    subproblem [default: 100]
  --normK=<val>                     norm of BlockOperator (K). Normally
                                    calculated and printed.
                                    Can then be reused to save time.
  --onlyNormK                       Calculate normK and then exit.
  --precond                         Use preconditioning
  --numSegsToCombine=<val>          Rebin all sinograms, with a given number of
                                    segments to combine. Increases speed.
  --numViewsToCombine=<val>         Rebin all sinograms, with a given number of
                                    views to combine. Increases speed.
  --normaliseDataAndBlock           Normalise raw data and block operator by
                                    multiplying by 1./normK.
  --no_log=<str>                    Disable log file.
  --algorithm=<string>              Which algorithm to run [default: spdhg]
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

from functools import partial
from os import path
from glob import glob
from docopt import docopt
from sirf.Utilities import error, show_2D_array
import pylab
import sirf.Reg as reg
import sirf.STIR as pet
from ccpi.optimisation.algorithms import PDHG, SPDHG
from ccpi.optimisation.functions import \
    KullbackLeibler, BlockFunction, IndicatorBox
from ccpi.optimisation.operators import \
    CompositionOperator, BlockOperator, LinearOperator
from ccpi.plugins.regularisers import FGP_TV
from ccpi.filters import regularisers
import numpy as np

pet.AcquisitionData.set_storage_scheme('memory')

__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)


# Multiple files
trans_pattern = str(args['--trans']).replace('%', '*')
sino_pattern = str(args['--sino']).replace('%', '*')
attn_pattern = str(args['--attn']).replace('%', '*')
rand_pattern = str(args['--rand']).replace('%', '*')
num_iters = int(args['--iter'])
regularisation = args['--reg']
trans_type = args['--trans_type']


if attn_pattern is None:
    attn_pattern = ""
if rand_pattern is None:
    rand_pattern = ""

# Norm
norm_file = args['--norm']
if norm_file:
    if not path.isfile(norm_file):
        raise error("Norm file not found: " + norm_file)

# Number of voxels
nxny = int(args['--nxny'])

# Output file
outp_prefix = args['--outp']

# Initial estimate
initial_estimate = args['--initial']

visualisations = True if args['--visualisations'] else False
nifti = True if args['--nifti'] else False
use_gpu = True if args['--gpu'] else False
descriptive_fname = True if args['--descriptive_fname'] else False
update_obj_fn_interval = int(args['--update_obj_fn_interval'])

# Verbosity
verbosity = int(args['--verbosity'])
pet.set_verbosity(verbosity)
if verbosity == 0:
    msg_red = pet.MessageRedirector(None, None, None)

# Save interval
save_interval = int(args['--save_interval'])
save_interval = min(save_interval, num_iters)

# Convergence variables
r_alpha = float(args['--alpha'])
r_iters = float(args['--reg_iters'])
precond = True if args['--precond'] else False
# algorithm selection
algorithm = str(args['--algorithm'])


def get_resampler(image, ref=None, trans=None):
    """Return a NiftyResample object for the specified transform and image."""
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


def get_filenames():
    """Get filenames."""
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
    if len(trans_files) > 0 and num_ms != len(trans_files):
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

    return [num_ms, trans_files, sino_files, attn_files, rand_files]


def read_files(trans_files, sino_files, attn_files, rand_files):
    """Read files."""
    if not trans_files:
        trans = None
    else:
        if trans_type == "tm":
            trans = [reg.AffineTransformation(file) for file in trans_files]
        elif trans_type == "disp":
            trans = [reg.NiftiImageData3DDisplacement(file)
                     for file in trans_files]
        elif trans_type == "def":
            trans = [reg.NiftiImageData3DDeformation(file)
                     for file in trans_files]
        else:
            raise error("Unknown transformation type")

    sinos_raw = [pet.AcquisitionData(file) for file in sino_files]
    attns = [pet.ImageData(file) for file in attn_files]
    rands_raw = [pet.AcquisitionData(file) for file in rand_files]

    return [trans, sinos_raw, attns, rands_raw]


def pre_process_sinos(sinos_raw, num_ms):
    """Preprocess raw sinograms.

    Make positive if necessary and do any required rebinning."""
    # If empty (e.g., no randoms), return
    if not sinos_raw:
        return sinos_raw
    # Loop over all sinograms
    sinos = [0]*num_ms
    for ind in range(num_ms):
        # If any sinograms contain negative values
        # (shouldn't be the case), set them to 0
        sino_arr = sinos_raw[ind].as_array()
        if (sino_arr < 0).any():
            print("Input sinogram " + str(ind) +
                  " contains -ve elements. Setting to 0...")
            sinos[ind] = sinos_raw[ind].clone()
            sino_arr[sino_arr < 0] = 0
            sinos[ind].fill(sino_arr)
        else:
            sinos[ind] = sinos_raw[ind]
        # If rebinning is desired
        segs_to_combine = 1
        if args['--numSegsToCombine']:
            segs_to_combine = int(args['--numSegsToCombine'])
        views_to_combine = 1
        if args['--numViewsToCombine']:
            views_to_combine = int(args['--numViewsToCombine'])
        if segs_to_combine * views_to_combine > 1:
            sinos[ind] = sinos[ind].rebin(segs_to_combine, views_to_combine)
            # only print first time
            if ind == 0:
                print(f"Rebinned sino dimensions: {sinos[ind].dimensions()}")

    return sinos


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
            image.fill(0.0)

    return image


def resample_attn_images(num_ms, attns, trans):
    """Resample attenuation images if necessary."""
    resampled_attns = None
    if trans is None:
        resampled_attns = attns
    else:
        if len(attns) > 0:
            resampled_attns = [0]*num_ms
            # if using GPU, dimensions of attn and recon images have to match
            ref = image if use_gpu else None
            for i in range(num_ms):
                # if we only have 1 attn image, then we need to resample into
                # space of each gate. However, if we have num_ms attn images,
                # then assume they are already in the correct position, so use
                # None as transformation.
                tran = trans[i] if len(attns) == 1 else None
                # If only 1 attn image, then resample that. If we have num_ms
                # attn images, then use each attn image of each frame.
                attn = attns[0] if len(attns) == 1 else attns[i]
                resam = get_resampler(attn, ref=ref, trans=tran)
                resampled_attns[i] = resam.forward(attn)
    return resampled_attns


def set_up_acq_models(num_ms, sinos, rands, resampled_attns, image):
    """Set up acquisition models."""
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


def create_kls(sinos, etas, scale_factor=None):
    """Create list of KullbackLeiblers from list of sinos and etas."""
    if scale_factor:
        scaled_sinos = [sino.clone()*scale_factor for sino in sinos]
        scaled_etas = [eta.clone()*scale_factor for eta in etas]
    else:
        scaled_sinos = sinos
        scaled_etas = etas
    kl = [KullbackLeibler(b=sino, eta=eta)
          for sino, eta in zip(sinos, etas)]
    return kl


def norm(self):
    """Norm implementation."""
    return LinearOperator.PowerMethod(self, 10)[0]


def set_up_reconstructor(acq_models, resamplers, algorithm, sinos, rands=None):
    """Set up reconstructor."""
    # Create composition operators containing linear
    # acquisition models and resamplers
    if resamplers is None:
        C = [am.get_linear_acquisition_model() for am in acq_models]
        # Need an implementation of the norm
        setattr(pet.AcquisitionModel, 'norm', norm)
    else:
        C = [CompositionOperator(
            am.get_linear_acquisition_model(), res, preallocate=True)
                for am, res in zip(*(acq_models, resamplers))]

    # We'll need an additive term (eta). If randoms are present, use them
    # Else, use a scaled down version of the sinogram
    etas = rands if rands else [sino * 0 + 1e-5 for sino in sinos]

    # NormK only needed for PDHG, not SPDHG
    if algorithm == 'pdhg':
        # Configure the PDHG algorithm
        if args['--normK'] and not args['--onlyNormK']:
            normK = float(args['--normK'])
        else:
            kl = create_kls(sinos, etas)
            f = BlockFunction(*kl)
            K = BlockOperator(*C)
            # Calculate normK
            print("Calculating norm of the block operator...")
            normK = K.norm(iterations=10)
            print("Norm of the BlockOperator ", normK)
            if args['--onlyNormK']:
                exit(0)
    else:
        normK = None

    # Optionally rescale sinograms and BlockOperator using normK
    scale_factor = 1./normK if args['--normaliseDataAndBlock'] else None
    kl = create_kls(sinos, etas, scale_factor)
    f = BlockFunction(*kl)
    K = BlockOperator(*C)
    if scale_factor:
        K *= scale_factor

    return [f, K, normK]


def get_nonzero_recip(data):
    """Get the reciprocal of a datacontainer.

    Voxels where input == 0
    will have their reciprocal set to 1 (instead of infinity)
    """
    inv_np = data.as_array()
    inv_np[inv_np == 0] = 1
    inv_np = 1./inv_np
    data.fill(inv_np)


def precond_proximal(self, x, tau, out=None):
    """Modify proximal method to work with preconditioned tau."""
    pars = {'algorithm': FGP_TV,
            'input': np.asarray(x.as_array()/tau.as_array(),
                                dtype=np.float32),
            'regularization_parameter': self.lambdaReg,
            'number_of_iterations': self.iterationsTV,
            'tolerance_constant': self.tolerance,
            'methodTV': self.methodTV,
            'nonneg': self.nonnegativity,
            'printingOut': self.printing}

    res, info = \
        regularisers.FGP_TV(pars['input'],
                            pars['regularization_parameter'],
                            pars['number_of_iterations'],
                            pars['tolerance_constant'],
                            pars['methodTV'],
                            pars['nonneg'],
                            self.device)
    if out is not None:
        out.fill(res)
    else:
        out = x.copy()
        out.fill(res)
    out *= tau
    return out


def get_tau_sigma(normK):
    """Get tau and sigma.

    If normK is None (because not required, then return
    sigma and tau as None, too."""
    if not normK:
        return [None, None]
    if precond:

        tau = K.adjoint(K.range_geometry().allocate(1))
        get_nonzero_recip(tau)

        tmp_sigma = K.direct(K.domain_geometry().allocate(1))
        sigma = 0.*tmp_sigma
        get_nonzero_recip(sigma[0])

        FGP_TV.proximal = precond_proximal
        print("Will run proximal with preconditioned tau...")

    # If not preconditioned
    else:
        if args['--sigma']:
            sigma = float(args['--sigma'])
        else:
            sigma = 1.0/normK
        # If we need to calculate default tau
        if args['--tau']:
            tau = float(args['--tau'])
        else:
            tau = 1.0/normK

    return [tau, sigma]


def set_up_regularisation():
    """Set up regularisation."""
    if regularisation == 'none':
        G = IndicatorBox(lower=0)
    elif regularisation == 'FGP_TV':
        r_iters = float(args['--reg_iters'])
        r_tolerance = 1e-7
        r_iso = 0
        r_nonneg = 1
        r_printing = 0
        device = 'gpu' if use_gpu else 'cpu'
        G = FGP_TV(r_alpha, r_iters, r_tolerance,
                   r_iso, r_nonneg, r_printing, device)
    else:
        raise error("Unknown regularisation")

    return G


def get_output_filename(attn_files, normK, sigma, tau, sino_files, resamplers):
    """Get output filename."""

    outp_file = outp_prefix
    if descriptive_fname:
        if len(attn_files) > 0:
            outp_file += "_wAC"
        if norm_file:
            outp_file += "_wNorm"
        if args['--rand']:
            outp_file += "_wRands"
        if use_gpu:
            outp_file += "_wGPU"
        if args['--normaliseDataAndBlock']:
            outp_file += '_wDataScale'
        else:
            outp_file += '_noDataScale'
        if algorithm == 'pdhg':
            if args['--normK']:
                outp_file += '_userNormK' + str(normK)
            else:
                outp_file += '_computedNormK' + str(normK)
            if not precond:
                outp_file += "_sigma" + str(sigma)
                outp_file += "_tau" + str(tau)
            else:
                outp_file += "_wPrecond"
        outp_file += "_Reg-" + regularisation
        if regularisation == 'FGP_TV':
            outp_file += "-alpha" + str(r_alpha)
            outp_file += "-riters" + str(r_iters)
        outp_file += '_' + algorithm
        outp_file += "_nGates" + str(len(sino_files))
        if resamplers is None:
            outp_file += "_noMotion"
    return outp_file


def get_algo(f, G, K, sigma, tau, outp_file):
    """Get the reconstruction algorithm."""
    if algorithm == 'pdhg':

        Algo = partial(PDHG, sigma=sigma, tau=tau)

    elif algorithm == 'spdhg':
        # let's define the subsets as the motion states
        num_subsets = len(K)
        # assign the probabilities implicit form
        prob = [1/num_subsets]*num_subsets
        # assign the probabilities explicit form
        # prob = [(num_subsets-1)*1/(2*num_subsets)] + [1/2]

        Algo = partial(SPDHG, prob=prob, sigma=None, tau=None)

    else:
        raise error("Unknown algorithm: " + algorithm)

    algo = Algo(f=f, g=G, operator=K,
                max_iteration=num_iters,
                update_objective_interval=update_obj_fn_interval,
                log_file=outp_file+".log",
                use_axpby=False)

    return algo


def get_save_callback_function(outp_file):
    """Get the save callback function."""
    def save_callback(save_interval, nifti, outp_file,
                      num_iters, iteration,
                      last_objective, x):
        """Save callback function."""
        completed_iterations = iteration + 1
        if completed_iterations % save_interval == 0 or \
                completed_iterations == num_iters:
            if not nifti:
                x.write("{}_iters_{}".format(outp_file, completed_iterations))
            else:
                reg.NiftiImageData(x).write(
                    "{}_iters_{}".format(outp_file, completed_iterations))

    psave_callback = partial(
        save_callback, save_interval, nifti, outp_file, num_iters)
    return psave_callback


def display_results(out_arr):
    """Display results if desired."""
    if visualisations:
        # show reconstructed image
        out_arr = algo.get_output().as_array()
        z = out_arr.shape[0]//2
        show_2D_array('Reconstructed image', out_arr[z, :, :])
        pylab.show()


def main():
    """Run main function."""
    ###########################################################################
    # Parse input files
    ###########################################################################

    [num_ms, trans_files, sino_files, attn_files, rand_files] = get_filenames()

    ###########################################################################
    # Read input
    ###########################################################################

    [trans, sinos_raw, attns, rands_raw] = \
        read_files(trans_files, sino_files, attn_files, rand_files)

    sinos = pre_process_sinos(sinos_raw, num_ms)
    rands = pre_process_sinos(rands_raw, num_ms)

    ###########################################################################
    # Initialise recon image
    ###########################################################################

    image = get_initial_estimate(sinos)

    ###########################################################################
    # Set up resamplers
    ###########################################################################

    if trans is None:
        resamplers = None
    else:
        resamplers = [get_resampler(image, trans=tran) for tran in trans]

    ###########################################################################
    # Resample attenuation images (if necessary)
    ###########################################################################

    resampled_attns = resample_attn_images(num_ms, attns, trans)

    ###########################################################################
    # Set up acquisition models
    ###########################################################################

    acq_models = set_up_acq_models(
        num_ms, sinos, rands, resampled_attns, image)

    ###########################################################################
    # Set up reconstructor
    ###########################################################################

    [f, K, normK] = set_up_reconstructor(
        acq_models, resamplers, algorithm, sinos, rands)

    # Get tau and sigma (scalars if no preconditioning, else they'll be arrays)
    [tau, sigma] = get_tau_sigma(normK)

    # Set up regularisation
    G = set_up_regularisation()

    # Get output filename
    outp_file = get_output_filename(
        attn_files, normK, sigma, tau, sino_files, resamplers)

    # Get algorithm
    algo = get_algo(f, G, K, sigma, tau, outp_file)

    # Create save call back function
    save_callback = get_save_callback_function(outp_file)

    # Run the reconstruction
    algo.run(num_iters, verbose=True, very_verbose=True,
             callback=save_callback)

    # Display results
    display_results(algo.get_output().as_array())


if __name__ == "__main__":
    main()
