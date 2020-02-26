# SIRF imports 
import sirf.Gadgetron as pMR
import sirf.STIR as pet
import sirf.Reg as reg
from sirf.Utilities import examples_data_path

# ccpi CIL imports
from ccpi.utilities.display import *
from ccpi.optimisation.algorithms import PDHG
from ccpi.optimisation.functions import LeastSquares, KullbackLeibler, BlockFunction, IndicatorBox
from ccpi.optimisation.operators import CompositionOperator, BlockOperator, LinearOperator
from ccpi.plugins.regularisers import FGP_TV

# import further modules
import os
import glob
import numpy as np

import matplotlib.pyplot as plt

pet.AcquisitionData.set_storage_scheme('memory')

#%%
# '''
# I've created 2D images of:  uMap, FDG, T1 and T2. They all have dimensions (1,150,150). 

# I've placed all of these into four motion states and resampled to give: 
# FDG_mf0.nii 
# uMap_mf0.nii 
# T1_mf0.nii 
# T2_mf0.nii 
# '''

def get_resampler_from_tm(tm, image):
    '''returns a NiftyResample object for the specified transform matrix and image'''

    mat = tm.as_array()

    resampler = reg.NiftyResample()
    resampler.set_reference_image(image)
    resampler.set_floating_image(image)
    resampler.add_transformation(tm)
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_linear()
    
    return resampler

def get_acquisition_model(uMap, templ_sino):
    '''create an acquisition model with uMap'''

    #%% create acquisition model
    am = pet.AcquisitionModelUsingRayTracingMatrix()
    am.set_num_tangential_LORs(5)

    # Set up sensitivity due to attenuation
    asm_attn = pet.AcquisitionSensitivityModel(uMap, am)
    asm_attn.set_up(templ_sino)
    bin_eff = pet.AcquisitionData(templ_sino)
    bin_eff.fill(1.0)
    print('applying attenuation (please wait, may take a while)...')
    asm_attn.unnormalise(bin_eff)
    asm_attn = pet.AcquisitionSensitivityModel(bin_eff)

    am.set_acquisition_sensitivity(asm_attn)

    am.set_up(templ_sino,uMap)
    return am

def convert_nifti_to_stir_ImageData(image, templ_sino, dim):
    # pet_target_image_template_reg = reg.NiftiImageData(os.path.basename(FDG_files[0]))
    out = pet.ImageData(templ_sino)
    #dim=(1,150,150)
    voxel_size=out.voxel_sizes()
    out.initialise(dim,voxel_size)
    out.fill(image.as_array())
    return out

data_path = os.path.join( os.path.abspath('/home/ofn77899') ,
                          'brainweb', 'brainweb_single_slice')

# number of motion states
n_ms = 4

FDG_files  = sorted( glob.glob(os.path.join( data_path, 'FDG_mf*.nii') ) )
uMap_files = sorted( glob.glob(os.path.join( data_path, 'uMap_mf*.nii') ) )
T1_files   = sorted( glob.glob(os.path.join( data_path, 'T1_mf*.nii') ) )
T2_files   = sorted( glob.glob(os.path.join( data_path, 'T2_mf*.nii') ) )
transform_matrices_files  = sorted( glob.glob(os.path.join( data_path, 'fwd_tm*.txt') ) )


# Acquisition Data
noisy_sinogram_files = sorted( glob.glob(os.path.join (data_path, 'sino_FDG_noisy_mf*.hs' )) )

# 

rotated_sinos = []
resamplers = []
ams = []

# needs to go to the data directory
os.chdir(data_path)
dim = (1, 150, 150)
# We'll need a template sinogram
mmr_template_sino_single_slice = os.path.join( data_path ,'mmr_single_slice_template.hs')
templ_sino = pet.AcquisitionData(mmr_template_sino_single_slice)
pet_target_image_template_reg = reg.NiftiImageData(os.path.basename(FDG_files[0]))
pet_target_image_template = convert_nifti_to_stir_ImageData(pet_target_image_template_reg, templ_sino, dim)

for ms in range(n_ms):
    # load the noisy Sinogram
    rotated_sinos.append( 
       pet.AcquisitionData(os.path.basename(noisy_sinogram_files[ms]))
    )
    # create the resamplers from the TransformationMatrix
    tm = reg.AffineTransformation(os.path.basename(transform_matrices_files[ms]))
    resamplers.append(

        get_resampler_from_tm(tm, pet_target_image_template)
 
    )
    # create the Acquisition Model with appropriate uMap
    # because of https://github.com/CCPPETMR/SIRF/issues/562 we 
    # pass from NumPy
    uMap = reg.NiftiImageData(os.path.basename(uMap_files[ms]))
    # im = pet.ImageData(templ_sino)
    # dim=(1,150,150)
    # voxel_size=im.voxel_sizes()
    # im.initialise(dim,voxel_size)
    # im.fill(uMap.as_array())
    im = convert_nifti_to_stir_ImageData(uMap, templ_sino, dim)

    ams.append(
        get_acquisition_model(im, templ_sino)
    )

# compose the resampler with the acquisition models
C = [ CompositionOperator(am, res) for am, res in zip (*(ams, resamplers)) ]

# the first AcquisitionModel doesn't need to be resampled we will try to reconstruct in 
# the orientation of the motion state 0.
# C[0] = ams[0]

# Configure the PDHG algorithm

kl = [ KullbackLeibler(b=rotated_sino, eta=(rotated_sino * 0 + 1e-5)) for rotated_sino in rotated_sinos ] 
f = BlockFunction(*kl)
K = BlockOperator(*C)
normK = K.norm(iterations=10)
#default values
sigma = 1/normK
tau = 1/normK 

print ("Norm of the BlockOperator ", normK)

    
# TV regularisation
#regularisation parameters for TV
# 
r_alpha = 5e-1
r_iterations = 500
r_tolerance = 1e-7
r_iso = 0
r_nonneg = 1
r_printing = 0

TV = FGP_TV(r_alpha, r_iterations, r_tolerance, r_iso,r_nonneg,r_printing,'cpu')

print ("current dir", os.getcwd())

G = IndicatorBox(lower=0)

pdhg = PDHG(f = f, g = G, operator = K, sigma = sigma, tau = tau, 
            max_iteration = 1000,
            update_objective_interval = 1)

pdhg.run(100, verbose=False)

solution = reg.NiftiImageData(os.path.basename(FDG_files[0]))
plotter2D([solution.as_array(), pdhg.get_output().as_array()[0]] )
# %%
