"""MCIR for PET

Usage:
  PET_MCIR [--help | options]

Options:
  -T <pattern>, --trans=<pattern>     transformation pattern (e.g., tm_ms*.txt). Enclose in quotations.
  -t <str>, --trans_type=<str>        transformation type (tm, disp, def) [default: tm]
  -S <pattern>, --sino=<pattern>      sinogram pattern (e.g., sino_ms*.hs). Enclose in quotations.
  -a <pattern>, --attn=<pattern>      attenuation pattern (e.g., attn_ms*.hv). Enclose in quotations.
  -R <pattern>, --rand=<pattern>      randoms pattern (e.g., rand_ms*.hs). Enclose in quotations.
  -n <norm>, --norm=<norm>            ECAT8 bin normalization file
  -i <int>, --iter=<int>              num iterations [default: 10]
  -r <string>, --reg=<string>         regularisation ("none","FGP_TV", ...) [default: none]
  -o <outp>, --outp=<outp>            output file prefix [default: recon]
  -d <nxny>, --nxny=<nxny>            image x and y dimensions as string '(nx,ny)'
                                      (no space after comma) [default: (127,127)]
  --visualisations                    show visualisations
  --nifti                             save output as nifti
  -v <int>, --verbosity=<int>         STIR verbosity [default: 0]
  -s <int>, --save_interval=<int>     save every x iterations [default: 10]
"""

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2020 University College London.
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
from ast import literal_eval
from glob import glob
from docopt import docopt
from sirf.Utilities import error, show_2D_array
import pylab
import sirf.Reg as reg
import sirf.STIR as pet
from ccpi.optimisation.algorithms import PDHG
from ccpi.optimisation.functions import KullbackLeibler, BlockFunction, IndicatorBox
from ccpi.optimisation.operators import CompositionOperator, BlockOperator
from ccpi.plugins.regularisers import FGP_TV

pet.AcquisitionData.set_storage_scheme('memory')

__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)


def file_exists(filename):
    """Check if file exists, optionally throw error if not"""
    return os.path.isfile(filename)


def check_file_exists(filename):
    """Check file exists, else throw error"""
    if not file_exists:
        raise error('File not found: %s' % filename)


# Multiple files
trans_pattern = str(args['--trans'])
sino_pattern = str(args['--sino'])
attn_pattern = str(args['--attn'])
rand_pattern = str(args['--rand'])
num_iters = int(args['--iter'])
regularisation = str(args['--reg'])
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
        raise error("Norm file not found: " + norm_file)

# Number of voxels
nxny = literal_eval(args['--nxny'])

# Output file
outp_prefix = str(args['--outp'])

if args['--visualisations']:
    visualisations = True
else:
    visualisations = False

if args['--nifti']:
    nifti = True
else:
    nifti = False

# Verbosity
pet.set_verbosity(int(args['--verbosity']))

# Verbosity
save_interval = int(args['--save_interval'])


def get_resampler_from_trans(trans, image):
    """returns a NiftyResample object for the specified transform and image"""
    resampler = reg.NiftyResample()
    resampler.set_reference_image(image)
    resampler.set_floating_image(image)
    resampler.add_transformation(trans)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
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

    ############################################################################################
    # Parse input files
    ############################################################################################

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
                             "#sinos = " + str(num_ms) + ", #trans = " + str(len(trans_files)))
    # If any rand, check num == num_ms
    if len(rand_files) > 0 and len(rand_files) != num_ms:
        raise AssertionError("#rand should match #sinos. "
                             "#sinos = " + str(num_ms) + ", #rand = " + str(len(rand_files)))

    # For attn, there should be 0, 1 or num_ms images
    if len(attn_files) != 0 and len(attn_files) != 1 and len(attn_files) != num_ms:
        raise AssertionError("#attn should be 0, 1 or #sinos")
    
    ############################################################################################
    # Read input
    ############################################################################################

    if trans_type == "tm":
        trans = [reg.AffineTransformation(file) for file in trans_files]
    elif trans_type == "disp":
        trans = [reg.NiftiImageData3DDisplacement(file) for file in trans_files]
    elif trans_type == "def":
        trans = [reg.NiftiImageData3DDeformation(file) for file in trans_files]
    else:
        raise error("Unknown transformation type")

    sinos_raw = [pet.AcquisitionData(file) for file in sino_files]
    attns = [pet.ImageData(file) for file in attn_files]
    rands = [pet.AcquisitionData(file) for file in rand_files]

    # If any sinograms contain negative values (shouldn't be the case), set them to 0
    sinos = [0]*num_ms
    for ind in range(num_ms):
        if (sinos_raw[ind].as_array() < 0).any():
            print("Input sinogram " + str(ind) + " contains -ve elements. Setting to 0...")
            sinos[ind] = sinos_raw[ind].clone()
            sino_arr = sinos[ind].as_array()
            sino_arr[sino_arr < 0] = 0
            sinos[ind].fill(sino_arr)
        else:
            sinos[ind] = sinos_raw[ind]

    ############################################################################################
    # Initialise recon image
    ############################################################################################

    image = sinos[0].create_uniform_image(1.0, nxny)

    ############################################################################################
    # Set up resamplers
    ############################################################################################

    resamplers = [get_resampler_from_trans(tran, image) for tran in trans]

    ############################################################################################
    # Set up acquisition models
    ############################################################################################

    print("Setting up acquisition models...")
    acq_models = num_ms * [pet.AcquisitionModelUsingRayTracingMatrix()]

    # If present, create ASM from ECAT8 normalisation data
    asm_norm = None
    if norm_file:
        asm_norm = pet.AcquisitionSensitivityModel(norm_file)

    # Loop over each motion state
    for ind in range(num_ms):
        # Create attn ASM if necessary
        asm_attn = None
        if len(attns) == num_ms:
            asm_attn = get_asm_attn(sinos[ind], attns[ind], acq_models[ind])
        elif len(attns) == 1:
            print("Resampling attn im " + str(ind) + " into required motion state...")
            resampler = get_resampler_from_trans(trans[ind], attns[0])
            resampled_attn = resampler.forward(attns[0])
            asm_attn = get_asm_attn(sinos[ind], resampled_attn, acq_models[ind])

        # Get ASM dependent on attn and/or norm
        asm = None
        if asm_norm and asm_attn:
            print("AcquisitionSensitivityModel contains norm and attenuation...")
            asm = pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
        elif asm_norm:
            print("AcquisitionSensitivityModel contains norm...")
            asm = asm_norm
        elif asm_attn:
            print("AcquisitionSensitivityModel contains attenuation...")
            asm = asm_attn
        if asm:
            print("Setting AcquisitionSensitivityModel...")
            acq_models[ind].set_acquisition_sensitivity(asm)

        if len(rands) > 0:
            acq_models[ind].set_background_term(rands[ind])

        # Set up
        acq_models[ind].set_up(sinos[ind], image)

    ############################################################################################
    # Set up reconstructor
    ############################################################################################

    print("Setting up reconstructor...")

    # Create composition operators containing acquisition models and resamplers
    C = [ CompositionOperator(am, res, preallocate=True) for am, res in zip (*(acq_models, resamplers)) ]

    # Configure the PDHG algorithm
    kl = [ KullbackLeibler(b=sino, eta=(sino * 0 + 1e-5)) for sino in sinos ] 
    f = BlockFunction(*kl)
    K = BlockOperator(*C)
    normK = K.norm(iterations=10)

    # normK = LinearOperator.PowerMethod(K, iterations=10)[0]
    # default values
    sigma = 1/normK
    tau = 1/normK 
    sigma = 0.001
    tau = 1/(sigma*normK**2)
    print("Norm of the BlockOperator ", normK)

    if regularisation == 'none':
        G = IndicatorBox(lower=0)
    elif regularisation == 'FGP_TV':
        r_alpha = 5e-1
        r_iterations = 100
        r_tolerance = 1e-7
        r_iso = 0
        r_nonneg = 1
        r_printing = 0
        G = FGP_TV(r_alpha, r_iterations, r_tolerance, r_iso, r_nonneg, r_printing, 'gpu')
    else:
        raise error("Unknown regularisation")

    pdhg = PDHG(f=f, g=G, operator=K, sigma=sigma, tau=tau,
                max_iteration=1000,
                update_objective_interval=1)

    # Get filename
    outp_file = outp_prefix
    if len(attn_files) > 0:
        outp_file += "_wAC"
    if norm_file:
        outp_file += "_wNorm"
    outp_file += "_Reg-" + regularisation
    outp_file += "_nGates" + str(len(sino_files))

    for i in range(1, num_iters+1, save_interval):
        pdhg.run(save_interval, verbose=True)
        out = pdhg.get_output()    
        if not nifti:
            out.write(outp_file + "_iters" + str(i))
        else:
            reg.NiftiImageData(out).write(outp_file + "_iters" + str(i))

    if visualisations:
        # show reconstructed image
        out_arr = out.as_array()
        z = out_arr.shape[0]//2
        show_2D_array('Reconstructed image', out.as_array()[z, :, :])
        pylab.show()


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
