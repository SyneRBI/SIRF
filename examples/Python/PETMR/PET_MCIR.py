"""MCIR for PET

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
  --sigma=<val>                     set PDHG sigma [default: 0.001]
  --tau=<val>                       set PDHG tau (default: 1/(sigma*normK**2)
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

import os
from glob import glob
from docopt import docopt
from sirf.Utilities import error, show_2D_array
import pylab
import sirf.Reg as reg
import sirf.STIR as pet
from ccpi.optimisation.algorithms import PDHG
from ccpi.optimisation.functions import KullbackLeibler, BlockFunction
from ccpi.optimisation.functions import IndicatorBox
from ccpi.optimisation.operators import CompositionOperator, BlockOperator
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
    if not os.path.isfile(norm_file):
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
precond = True if args['--precond'] else False


def get_resampler(image, ref=None, trans=None):
    """returns a NiftyResample object for the specified transform and image"""
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
    """Get attn ASM from sino, attn image and acq model"""
    asm_attn = pet.AcquisitionSensitivityModel(attn, acq_model)
    # temporary fix pending attenuation offset fix in STIR:
    # converting attenuation into 'bin efficiency'
    asm_attn.set_up(sino)
    bin_eff = pet.AcquisitionData(sino)
    bin_eff.fill(1.0)
    asm_attn.unnormalise(bin_eff)
    asm_attn = pet.AcquisitionSensitivityModel(bin_eff)
    return asm_attn


def main():

    ###########################################################################
    # Parse input files
    ###########################################################################

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

    ###########################################################################
    # Read input
    ###########################################################################

    if trans_type == "tm":
        trans = [reg.AffineTransformation(file) for file in trans_files]
    elif trans_type == "disp":
        trans = [reg.NiftiImageData3DDisplacement(file)
                 for file in trans_files]
    elif trans_type == "def":
        trans = [reg.NiftiImageData3DDeformation(file) for file in trans_files]
    else:
        raise error("Unknown transformation type")

    sinos_raw = [pet.AcquisitionData(file) for file in sino_files]
    attns = [pet.ImageData(file) for file in attn_files]
    rands = [pet.AcquisitionData(file) for file in rand_files]

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

    ###########################################################################
    # Initialise recon image
    ###########################################################################

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

    ###########################################################################
    # Set up resamplers
    ###########################################################################

    resamplers = [get_resampler(image, trans=tran) for tran in trans]

    ###########################################################################
    # Resample attenuation images (if necessary)
    ###########################################################################

    resampled_attns = None
    if len(attns) > 0:
        resampled_attns = [0]*num_ms
        # if using GPU, dimensions of attn and recon images have to match
        ref = image if use_gpu else None
        for i in range(len(attns)):
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

    ###########################################################################
    # Set up acquisition models
    ###########################################################################

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
            asm_attn = get_asm_attn(sinos[ind], resampled_attns[i],
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

    ###########################################################################
    # Set up reconstructor
    ###########################################################################

    print("Setting up reconstructor...")

    # Create composition operators containing acquisition models and resamplers
    C = [CompositionOperator(am, res, preallocate=True)
         for am, res in zip(*(acq_models, resamplers))]

    # Configure the PDHG algorithm
    if args['--normK'] and not args['--onlyNormK']:
        normK = float(args['--normK'])
    else:
        kl = [KullbackLeibler(b=sino, eta=(sino * 0 + 1e-5)) for sino in sinos]
        f = BlockFunction(*kl)
        K = BlockOperator(*C)
        # Calculate normK
        print("Calculating norm of the block operator...")
        normK = K.norm(iterations=10)
        print("Norm of the BlockOperator ", normK)
        if args['--onlyNormK']:
            exit(0)

    # Optionally rescale sinograms and BlockOperator using normK
    scale_factor = 1./normK if args['--normaliseDataAndBlock'] else 1.0
    kl = [KullbackLeibler(b=sino*scale_factor, eta=(sino * 0 + 1e-5))
          for sino in sinos]
    f = BlockFunction(*kl)
    K = BlockOperator(*C)*scale_factor

    # If preconditioned
    if precond:

        def get_nonzero_recip(data):
            """Get the reciprocal of a datacontainer. Voxels where input == 0
            will have their reciprocal set to 1 (instead of infinity)"""
            inv_np = data.as_array()
            inv_np[inv_np == 0] = 1
            inv_np = 1./inv_np
            data.fill(inv_np)

        tau = K.adjoint(K.range_geometry().allocate(1))
        get_nonzero_recip(tau)

        tmp_sigma = K.direct(K.domain_geometry().allocate(1))
        sigma = 0.*tmp_sigma
        get_nonzero_recip(sigma[0])

        def precond_proximal(self, x, tau, out=None):
            """Modify proximal method to work with preconditioned tau"""
            pars = {'algorithm': FGP_TV,
                    'input': np.asarray(x.as_array()/tau.as_array(),
                                        dtype=np.float32),
                    'regularization_parameter': self.lambdaReg,
                    'number_of_iterations': self.iterationsTV,
                    'tolerance_constant': self.tolerance,
                    'methodTV': self.methodTV,
                    'nonneg': self.nonnegativity,
                    'printingOut': self.printing}

            res, info = regularisers.FGP_TV(pars['input'],
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

        FGP_TV.proximal = precond_proximal
        print("Will run proximal with preconditioned tau...")

    # If not preconditioned
    else:
        sigma = float(args['--sigma'])
        # If we need to calculate default tau
        if args['--tau']:
            tau = float(args['--tau'])
        else:
            tau = 1/(sigma*normK**2)

    if regularisation == 'none':
        G = IndicatorBox(lower=0)
    elif regularisation == 'FGP_TV':
        r_iterations = float(args['--reg_iters'])
        r_tolerance = 1e-7
        r_iso = 0
        r_nonneg = 1
        r_printing = 0
        device = 'gpu' if use_gpu else 'cpu'
        G = FGP_TV(r_alpha, r_iterations, r_tolerance,
                   r_iso, r_nonneg, r_printing, device)
    else:
        raise error("Unknown regularisation")

    if precond:
        def PDHG_new_update(self):
            """Modify the PDHG update to allow preconditioning"""
            # save previous iteration
            self.x_old.fill(self.x)
            self.y_old.fill(self.y)

            # Gradient ascent for the dual variable
            self.operator.direct(self.xbar, out=self.y_tmp)
            self.y_tmp *= self.sigma
            self.y_tmp += self.y_old

            self.f.proximal_conjugate(self.y_tmp, self.sigma, out=self.y)

            # Gradient descent for the primal variable
            self.operator.adjoint(self.y, out=self.x_tmp)
            self.x_tmp *= -1*self.tau
            self.x_tmp += self.x_old

            self.g.proximal(self.x_tmp, self.tau, out=self.x)

            # Update
            self.x.subtract(self.x_old, out=self.xbar)
            self.xbar *= self.theta
            self.xbar += self.x

        PDHG.update = PDHG_new_update

    # Get filename
    outp_file = outp_prefix
    if descriptive_fname:
        if len(attn_files) > 0:
            outp_file += "_wAC"
        if norm_file:
            outp_file += "_wNorm"
        if use_gpu:
            outp_file += "_wGPU"
        outp_file += "_Reg-" + regularisation
        if regularisation == 'FGP_TV':
            outp_file += "-alpha" + str(r_alpha)
            outp_file += "-riters" + str(r_iterations)
        if args['--normK']:
            outp_file += '_userNormK' + str(normK)
        else:
            outp_file += '_calcNormK' + str(normK)
        if args['--normaliseDataAndBlock']:
            outp_file += '_wDataScale'
        else:
            outp_file += '_noDataScale'
        if not precond:
            outp_file += "_sigma" + str(sigma)
            outp_file += "_tau" + str(tau)
        else:
            outp_file += "_wPrecond"
        outp_file += "_nGates" + str(len(sino_files))
        if resamplers is None:
            outp_file += "_noMotion"

    pdhg = PDHG(f=f, g=G, operator=K, sigma=sigma, tau=tau,
                max_iteration=num_iters,
                update_objective_interval=update_obj_fn_interval,
                x_init=image, log_file=outp_file+".log")

    def callback_save(iteration, objective_value, solution):
        """Callback function to save images"""
        if (iteration+1) % save_interval == 0:
            out = solution if not nifti else reg.NiftiImageData(solution)
            out.write(outp_file + "_iters" + str(iteration+1))

    pdhg.run(iterations=num_iters, callback=callback_save,
             verbose=True, very_verbose=True)

    if visualisations:
        # show reconstructed image
        out = pdhg.get_output()
        out_arr = out.as_array()
        z = out_arr.shape[0]//2
        show_2D_array('Reconstructed image', out.as_array()[z, :, :])
        pylab.show()


# if anything goes wrong, an exception will be thrown
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
    exit(0)
except error as err:
    # display error information
    print('%s' % err.value)
    exit(1)
