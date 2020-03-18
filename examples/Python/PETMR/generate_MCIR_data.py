"""Generate data for MCIR

Usage:
  generate_MCIR_data [--help | options]

Options:
  -p <path>, --path=<path>     path to output data [default: brainweb]
  -P <PET>, --PET=<PET>        template PET raw data (default: data/examples/PET/mMR/mMR_template_span11.hs)
  -M <MR>, --MR=<MR>           template MR raw data (default: data/examples/MR/grappa2_1rep.h5)
  -n <int>, --num_ms=<int>     number of motion states [default: 4]
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
import sirf.Gadgetron as mr
import sirf.STIR as pet
import sirf.Reg as reg
from sirf.Utilities import examples_data_path
import brainweb
from math import cos, sin, radians
import numpy as np
from random import randint
from sirf.Utilities import error
from docopt import docopt

__version__ = '0.1.0'
args = docopt(__doc__, version=__version__)

# path
data_path = args['--path']

# PET sino
if args['--PET']:
    template_PET_raw_path = args['--PET']
else:
    template_PET_raw_path = os.path.join(examples_data_path('PET'), 'mMR/mMR_template_single_slice.hs')

# MR sino
if args['--MR']:
    template_MR_raw_path = args['--MR']
else:
    template_MR_raw_path = os.path.join(examples_data_path('MR'), 'grappa2_1rep.h5')

# Num motion states
num_ms = int(args['--num_ms'])


def get_resampler_from_tm(tm, image):
    """returns a NiftyResample object for the specified transform matrix and image"""
    resampler = reg.NiftyResample()
    resampler.set_reference_image(image)
    resampler.set_floating_image(image)
    resampler.add_transformation(tm)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    
    return resampler


def simple_mr_recon(input_data):
    """Simple MR recon from input data"""
    recon = mr.CartesianGRAPPAReconstructor()
    recon.set_input(input_data)
    recon.compute_gfactors(False)
    recon.process()
    return recon.get_output()


def download_data():
    """Download brainweb data"""
    fname, url = sorted(brainweb.utils.LINKS.items())[0]
    files = brainweb.get_file(fname, url, ".")
    data = brainweb.load_file(fname)

    brainweb.seed(1337)

    vol = brainweb.get_mmr_fromfile(
        fname,
        petNoise=1, t1Noise=0.75, t2Noise=0.75,
        petSigma=1, t1Sigma=1, t2Sigma=1)

    FDG = vol['PET']
    uMap = vol['uMap']
    T1 = vol['T1']
    return [FDG, uMap, T1]


def crop_brainweb(template_im, vol, z_slice, xy_min, xy_max):
    """Crop from (127,344,344) to (1,nxy,nxy)"""
    im = template_im.clone()
    vol = vol[z_slice,xy_min:xy_max,xy_min:xy_max]
    im.fill(vol)
    return im


def get_and_save_tm(i):
    """Get and save affine transformation matrix"""
    if i == 0:
        [r, t_x, t_y] = [0., 0., 0.]
    elif i == 1: 
        [r, t_x, t_y] = [10., -10., 0.]
    elif i == 2: 
        [r, t_x, t_y] = [20., -5., 5.]
    elif i == 3: 
        [r, t_x, t_y] = [-10., 10., 5.]
    else:
        [r, t_x, t_y] = [randint(-20,20), randint(-20,20), randint(-20,20)]
        
    r = radians(r)
    tm = reg.AffineTransformation(np.array(
        [[ cos(r), sin(r), 0, t_x],
         [-sin(r), cos(r), 0, t_y],
         [      0,      0, 1, 0  ],
         [      0,      0, 0, 1  ]]))
    tm.write('fwd_tm_ms_' + str(i))
    return tm


def get_acquisition_model(uMap, templ_sino):
    """Create acquisition model"""
    am = pet.AcquisitionModelUsingRayTracingMatrix()
    am.set_num_tangential_LORs(5)

    # Set up sensitivity due to attenuation
    asm_attn = pet.AcquisitionSensitivityModel(uMap, am)
    asm_attn.set_up(templ_sino)
    bin_eff = pet.AcquisitionData(templ_sino)
    bin_eff.fill(1.0)
    asm_attn.unnormalise(bin_eff)
    asm_attn = pet.AcquisitionSensitivityModel(bin_eff)
    am.set_acquisition_sensitivity(asm_attn)
    am.set_up(templ_sino,uMap)
    return am


def add_noise(fraction_of_counts, sinogram):
    """Add noise to sinogram"""
    sino_arr = sinogram.as_array()
    minmax = (sino_arr.min(), sino_arr.max())
    if fraction_of_counts > 0 and fraction_of_counts < 1:
        fraction_of_counts = fraction_of_counts * (minmax[1] - minmax[0])
    else:
        raise AssertionError("Fraction of counts should be > 0 and < 1")
       
    sino_arr = fraction_of_counts * ((sino_arr-minmax[0]) / (minmax[1]-minmax[0]))
    noisy_sino = sinogram * 0.
    noisy_sino.fill( np.random.poisson(sino_arr) )
    return noisy_sino


def main():

    # Make output folder if necessary
    if not os.path.isdir(data_path):
        os.makedirs(data_path)
    os.chdir(data_path)

    # Download the data
    print("downloading brainweb data...")
    [FDG_arr, uMap_arr, T1_arr] = download_data()

    # Get template PET image from template raw
    template_PET_raw = pet.AcquisitionData(template_PET_raw_path)
    template_PET_im = pet.ImageData(template_PET_raw)

    # Get template MR image from template raw
    template_MR_raw = mr.AcquisitionData(template_MR_raw_path)
    template_MR_raw.sort_by_time()
    template_MR_raw = mr.preprocess_acquisition_data(template_MR_raw)
    template_MR_im = simple_mr_recon(template_MR_raw)

    # Number voxels in (x,y) directions - nxy (dictated by MR image)
    nxy = template_MR_im.get_geometrical_info().get_size()[0]
    if nxy != template_MR_im.get_geometrical_info().get_size()[1]:
        raise AssertionError("Expected square image in (x,y) direction")
    if template_MR_im.get_geometrical_info().get_size()[2] > 1:
        raise AssertionError("Only currently designed for 2D image")

    # Create PET image
    dim=(1,nxy,nxy)
    size = FDG_arr.shape
    z_slice = size[0] // 2
    xy_min = (size[1] - nxy) // 2
    xy_max = xy_min + nxy
    voxel_size=template_PET_im.voxel_sizes()
    template_PET_im.initialise(dim,voxel_size)

    # Reorient template MR image with template PET image such that it's compatible with both
    template_MR_im.reorient(template_PET_im.get_geometrical_info())

    ############################################################################################
    # Crop brainweb image to right size
    ############################################################################################

    # Convert brainweb's (127,344,344) to desired size
    print("Cropping brainweb images to size...")

    [FDG, uMap, T1] = [crop_brainweb(template_MR_im, im_arr, z_slice, xy_min, xy_max) \
        for im_arr in [FDG_arr, uMap_arr, T1_arr]]

    ############################################################################################
    # Apply motion
    ############################################################################################

    print("Resampling images to different motion states...")
    FDGs = [0]*num_ms
    uMaps = [0]*num_ms
    T1s = [0]*num_ms
    for ind in range(num_ms):
        # Get TM for given motion state
        tm = get_and_save_tm(ind)
        # Get resampler
        res = get_resampler_from_tm(tm, template_MR_im)
        # Resample
        for im, modality in zip([FDG, uMap, T1], ['FDG', 'uMap', 'T1']):
            resampled = res.forward(im)
            if modality == 'FDG':
                FDGs[ind] = resampled
            elif modality == 'uMap':
                uMaps[ind] = resampled
            elif modality == 'T1':
                T1s[ind] = resampled
            else:
                raise AssertionError("Unknown modality")
            reg.NiftiImageData(resampled).write(modality + '_ms' + str(ind))

    ############################################################################################
    # MR: create k-space data for motion states
    ############################################################################################

    # Create coil sensitivity data
    print("Calculating coil sensitivity map...")
    csm = mr.CoilSensitivityData()
    csm.smoothness = 500
    csm.calculate(template_MR_raw)

    # Create interleaved sampling
    print("Creating raw k-space data for MR motion states...")
    mvec = []
    for ind in range(num_ms):
        mvec.append(np.arange(ind, template_MR_raw.number(), num_ms))

    # Go through motion states and create k-space
    for ind in range(num_ms):

        acq_ms = template_MR_raw.new_acquisition_data(empty=True)

        # Set first two (??) acquisition
        acq_ms.append_acquisition(template_MR_raw.acquisition(0))
        acq_ms.append_acquisition(template_MR_raw.acquisition(1))

        # Add motion resolved data
        for jnd in range(len(mvec[ind])):

            if mvec[ind][jnd] < template_MR_raw.number() - 1 and mvec[ind][jnd] > 1:  # Ensure first and last are not added twice
                cacq = template_MR_raw.acquisition(mvec[ind][jnd])
                acq_ms.append_acquisition(cacq)

        # Set last acquisition
        acq_ms.append_acquisition(template_MR_raw.acquisition(template_MR_raw.number() - 1))

        # Create acquisition model
        AcqMod = mr.AcquisitionModel(acq_ms, T1s[ind])
        AcqMod.set_coil_sensitivity_maps(csm)

        # Forward project!
        acq_ms_sim = AcqMod.forward(T1s[ind])

        # Save
        print("writing: " + 'raw_T1_ms' + str(ind) + '.h5')
        acq_ms_sim.write('raw_T1_ms' + str(ind) + '.h5')

    ############################################################################################
    # PET: create sinograms
    ############################################################################################

    print("Creating singorams for PET motion states...")
    stir_uMap = template_PET_im.clone()
    stir_FDG = template_PET_im.clone()
    for ind in range(num_ms):
        stir_uMap.fill(uMaps[ind].as_array())
        stir_FDG.fill(FDGs[ind].as_array())
        am = get_acquisition_model(stir_uMap, template_PET_raw)
        FDG_sino = am.forward(stir_FDG)
        FDG_sino = add_noise(0.25, FDG_sino)
        FDG_sino.write('raw_FDG_ms' + str(ind))


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('done')
except error as err:
    # display error information
    print('%s' % err.value)
