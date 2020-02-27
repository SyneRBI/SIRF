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

import subprocess
import fileinput

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



def read_2D_STIR_nii(filename):
    temp_filename = filename + "_temp.hv"
    a = pet.ImageData(filename)
    a.write(temp_filename)
    text_to_search = "scaling factor (mm/pixel) [3] := 1"
    text_to_replace = "scaling factor (mm/pixel) [3] := 2.03125"
    dfile = fileinput.FileInput(temp_filename, inplace=True, backup='.bak')
    for line in dfile:
        print(line.replace(text_to_search, text_to_replace))
    dfile.close()
    im = pet.ImageData(temp_filename)
    return im

def add_noise(counts, sinogram):
    sino_arr = sinogram.as_array()
    minmax = (sino_arr.min(), sino_arr.max())
    if counts > 0 and counts <= 1:
        counts = counts * (minmax[1] - minmax[0])
    elif isinstance (counts, int):
        pass
       
    sino_arr = counts * ((sino_arr -minmax[0]) / (minmax[1]-minmax[0]))
    noisy_counts = sinogram * 0.
    noisy_counts.fill( np.random.poisson(sino_arr) )
    
    return noisy_counts



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
noiseless_sinogram_files = sorted( glob.glob(os.path.join (data_path, 'sino_FDG_mf*.hs' )) )
# 
# needs to go to the data directory
os.chdir(data_path)


#### 
## Shows the sinograms
if True:
    noisysino = pet.AcquisitionData(os.path.basename(noisy_sinogram_files[0])).as_array()
    sino = pet.AcquisitionData(os.path.basename(noiseless_sinogram_files[0])).as_array()
    sino_positive = pet.AcquisitionData(os.path.basename(noiseless_sinogram_files[0])).maximum(0).as_array()
    noisy_positive = add_noise(0.25, pet.AcquisitionData(os.path.basename(noiseless_sinogram_files[0])).maximum(0))
    ns_minmax = noisysino.min(), noisysino.max()
    s_minmax = sino.min(), sino.max()
    sp_minmax = sino_positive.min(), sino_positive.max()
    plotter2D( [sino[0][0], noisysino[0][0] , sino_positive[0][0], noisy_positive.as_array()[0][0]], 
    titles=['noiseless sino {}'.format(s_minmax), 'noisy sino {}'.format(ns_minmax), 
            'noiseless sino positive{}'.format(sp_minmax),
            'noisy positive 25%'] )


rotated_sinos = []
resamplers = []
ams = []


dim = (1, 150, 150)
# We'll need a template sinogram
mmr_template_sino_single_slice = os.path.join( data_path ,'mmr_single_slice_template.hs')
templ_sino = pet.AcquisitionData(mmr_template_sino_single_slice)

pet_target_image_template = read_2D_STIR_nii(os.path.basename(FDG_files[0]))

for ms in range(n_ms):
    # load the noisy Sinogram
    sino = pet.AcquisitionData(os.path.basename(noiseless_sinogram_files[ms])).maximum(0)
    
    rotated_sinos.append( 
        add_noise(0.25, sino)
    )
    # create the resamplers from the TransformationMatrix
    tm = reg.AffineTransformation(os.path.basename(transform_matrices_files[ms]))
    resamplers.append(

        get_resampler_from_tm(tm, pet_target_image_template)
 
    )
    # create the Acquisition Model with appropriate uMap
    uMap = read_2D_STIR_nii(os.path.basename(uMap_files[ms]))

    ams.append(
        get_acquisition_model(uMap, templ_sino)
    )

if True:
    # plot the solutions to see that with the resampler we get the same orientation
    implot = []
    for ms in range(n_ms):
        img = read_2D_STIR_nii(os.path.basename(FDG_files[ms]))
        implot.append(img.as_array()[0])
        implot.append(resamplers[ms].adjoint(img).as_array()[0])

    plotter2D( implot , titles=['MS0', 'Resampled to MS0', 
    'MS1', 'Resampled to MS0',
    'MS2', 'Resampled to MS0',
    'MS3', 'Resampled to MS0'])

# compose the resampler with the acquisition models
C = [ CompositionOperator(am, res) for am, res in zip (*(ams, resamplers)) ]
# C = [ am for am in ams ]
# C = [ ams[0] for i in range(len(ams))]
print ("number of motion states", len(resamplers))


# Configure the PDHG algorithm

kl = [ KullbackLeibler(b=rotated_sino, eta=(rotated_sino * 0 + 1e-5)) for rotated_sino in rotated_sinos ] 
f = BlockFunction(*kl)
K = BlockOperator(*C)

#f = kl[0]
#K = ams[0]

normK = K.norm(iterations=10)
#normK = LinearOperator.PowerMethod(K, iterations=10)[0]
#default values
sigma = 1/normK
tau = 1/normK 
sigma = 0.001
tau = 1/(sigma*normK**2)
print ("Norm of the BlockOperator ", normK)

    
# TV regularisation
#regularisation parameters for TV
# 
r_alpha = 5e-2
r_iterations = 100
r_tolerance = 1e-7
r_iso = 0
r_nonneg = 1
r_printing = 0

TV = FGP_TV(r_alpha, r_iterations, r_tolerance, r_iso,r_nonneg,r_printing,'gpu')

print ("current dir", os.getcwd())

G = IndicatorBox(lower=0)
# G = TV
pdhg = PDHG(f = f, g = G, operator = K, sigma = sigma, tau = tau, 
            max_iteration = 1000,
            update_objective_interval = 1)

pdhg.run(100, verbose=False)
#img = convert_nifti_to_stir_ImageData(reg.NiftiImageData(os.path.basename(FDG_files[0])), templ_sino, dim)
img = read_2D_STIR_nii(os.path.basename(FDG_files[0]))
solution = resamplers[0].direct(img)


# res = reg.NiftyResample()
# res.set_reference_image(pdhg.get_output())
# res.set_floating_image(solution)
# res.set_interpolation_type_to_linear()

# solution2 = res.direct(solution)

plotter2D([solution.as_array()[0], pdhg.get_output().as_array()[0]],
          titles = ['Ground Truth (MS0)' , 'PDHG output with TV'] )
# %%
