# -*- coding: utf-8 -*-
###
# Demonstration of advanced PET reconstruction with CCP PET-MR Software
# and randomized algorithms
#
# This demonstration shows how to use the stochastic primal-dual hybrid
# gradient algorithm (SPDHG) for regularized PET reconstruction.
#
#    [CERS2018] A. Chambolle, M. J. Ehrhardt, P. Richtarik and C.-B. Schoenlieb,
#    *Stochastic Primal-Dual Hybrid Gradient Algorithm with Arbitrary Sampling
#    and Imaging Applications*. SIAM Journal on Optimization, 28(4), 2783–2808
#    (2018) http://doi.org/10.1007/s10851-010-0251-1
#
#    [E+2017] M. J. Ehrhardt, P. J. Markiewicz, P. Richtarik, J. Schott,
#    A. Chambolle and C.-B. Schoenlieb, *Faster PET reconstruction with a
#    stochastic primal-dual hybrid gradient method*. Wavelets and Sparsity XVII,
#    58 (2017) http://doi.org/10.1117/12.2272946.
#
#    [EMS2018] M. J. Ehrhardt, P. J. Markiewicz and C.-B. Schoenlieb, *Faster
#    PET Reconstruction with Non-Smooth Priors by Randomization and
#    Preconditioning*. (2018) ArXiv: http://arxiv.org/abs/1808.07150
#
# This demo is a 'script', i.e. intended to be run step by step in a
# Python IDE such as spyder. It is organised in 'cells'. spyder displays these
# cells nicely and allows you to run each cell on its own.
#
# First version: 8th of September 2018
# Author: Matthias J Ehrhardt, Edoardo Pasca
#
## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2018 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2018 University College London.
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

# Initial imports etc
import os
import shutil
import numpy
import matplotlib.pyplot as plt
import pSTIR as pet
from ccpi.optimisation.spdhg import spdhg
from ccpi.optimisation.spdhg import KullbackLeibler
from ccpi.optimisation.funcs import ZeroFun
from ccpi.plugins.regularisers import FGP_TV, TGV, LLT_ROF, Diff4th
# from ccpi.optimisation.funcs import IndicatorBox
from spdhgutils import PowerMethodNonsquare
from spdhgutils import SubsetOperator
from spdhgutils import cilPluginToSIRFFactory

# plotting settings
# plt.ion() # interactive 'on' such that plots appear during loops
# some handy function definitions

def imshow3(image, sliceno=None, **kwargs):
    im = image.as_array()

    if sliceno is None:
        sliceno = im.shape[0]/2

    imshow(im[sliceno, :, :], **kwargs)

def imshow(image, newfig=True, limits=None, title=''):
    """Usage: imshow(image, [min,max], title)"""
    if newfig:
        plt.figure()

    plt.title(title)

    bitmap = plt.imshow(image)
    if limits is None:
        limits = [image.min(), image.max()]

    plt.clim(limits[0], limits[1])
    plt.colorbar(shrink=.6)
    plt.axis('off')
    return bitmap

# go to directory with input files
# adapt this path to your situation (or start everything in the relevant directory)
os.chdir(pet.petmr_data_path('pet'))
# copy files to working folder and change directory to where the output files are
EXAMPLE = 'nema'

if EXAMPLE == 'thorax':

    shutil.rmtree('working_folder/thorax_single_slice', True)
    shutil.copytree('thorax_single_slice', 'working_folder/thorax_single_slice')
    os.chdir('working_folder/thorax_single_slice')

elif EXAMPLE == 'brain':

    shutil.rmtree('working_folder/brain', True)
    shutil.copytree('brain', 'working_folder/brain')
    os.chdir('working_folder/brain')
elif EXAMPLE == 'nema':
    os.chdir('mMR')
    noisy_data = pet.AcquisitionData('sino_f1g1d0b0.hs')
else:
    raise ValueError('unexpected example')

if EXAMPLE in ['brain', 'thorax']:
    # We will first create some simulated data from ground-truth images

    # Read in images
    image = pet.ImageData('emission.hv')
    image_array = image.as_array()*.05
    image.fill(image_array)
    mu_map = pet.ImageData('attenuation.hv')
    mu_map_array = mu_map.as_array()
    # bitmap display of images
    sliceno = image_array.shape[0]/2
    plt.figure()
    imshow(image_array[sliceno, :, :, ], title='emission image')

    # save max for future displays
    cmax = image_array.max()*.6

    # create acquisition model
    am = pet.AcquisitionModelUsingRayTracingMatrix()
    am.set_num_tangential_LORs(5)
    templ = pet.AcquisitionData('template_sinogram.hs')
    am.set_up(templ, image)
    # simulate some data using forward projection
    acquired_data = am.forward(image)
    acquisition_array = acquired_data.as_array()

    # Display bitmaps of a middle sinogram
    plt.figure()
    imshow(acquisition_array[0, :, :, ], title='Forward projection')

    # close all plots
    plt.close('all')

    data = acquired_data
    background = data.copy()
    background.fill(5)

    array = data.as_array() + background.as_array()
    noisy_array = numpy.random.poisson(array.astype('float64'))
    max_counts = noisy_array.max()
    print(' Maximum counts in the data: {}'.format(noisy_array.max()))
    noisy_data = data.clone()
    noisy_data.fill(noisy_array)
elif EXAMPLE == 'nema':
    norm_file = 'norm.n.hdr'
    attn_file = 'mu_map.hv'
    sino_file = 'sino_f1g1d0b0.hs'
    rand_file = 'sino_randoms_f1g1d0b0.hs'
    noisy_data = pet.AcquisitionData(sino_file)
    randoms = pet.AcquisitionData(rand_file)
    asm_norm = pet.AcquisitionSensitivityModel(norm_file)

    am = pet.AcquisitionModelUsingRayTracingMatrix()
    am.set_num_tangential_LORs(10)
    # add it to the acquisition model
    am.set_acquisition_sensitivity(asm_norm)
    # Attenuation image
    attn_image = pet.ImageData(attn_file)
    asm_attn = pet.AcquisitionSensitivityModel(attn_image, am)
    # converting attenuation into attenuation factors (see previous exercise)
    asm_attn.set_up(noisy_data)
    attn_factors = pet.AcquisitionData(noisy_data)
    attn_factors.fill(1.0)
    print('applying attenuation (please wait, may take a while)...')
    asm_attn.unnormalise(attn_factors)
    print('done')
    # use these in the final attenuation model
    asm_attn = pet.AcquisitionSensitivityModel(attn_factors)
    # chain attenuation and normalisation
    asm = pet.AcquisitionSensitivityModel(asm_norm, asm_attn)
    # update the acquisition model etc
    am.set_acquisition_sensitivity(asm)

    am.set_background_term(randoms)
    # create initial image estimate of dimensions and voxel sizes
    # compatible with the scanner geometry (included in the AcquisitionData
    # object acq_data) and initialize each voxel to 1.0
    nxny = (127, 127)
    initial_image = noisy_data.create_uniform_image(1.0, nxny)
    am.set_up(noisy_data, initial_image)
    image = initial_image
    background = randoms


g_noreg = ZeroFun()

#%%

g_reg = cilPluginToSIRFFactory.getInstance(FGP_TV, 
                                           lambdaReg=.3,
                                           iterationsTV=1000,
                                           tolerance=1e-5,
                                           methodTV=0,
                                           nonnegativity=1,
                                           printing=0,
                                           device='cpu')
#%%
'''
g_reg = cilPluginToSIRFFactory.getInstance(LLT_ROF, 
                                           regularisation_parameterROF=0.04,
                                           regularisation_parameterLLT=0.01,
                                           iterations=500,
                                           time_marching_parameter=0.00002,
                                           device='cpu')

g_reg = cilPluginToSIRFFactory.getInstance(TGV, 
                                           regularisation_parameter=0.0005,
                                           alpha1=1,
                                           alpha0=0.7,
                                           iterations=250,
                                           LipshitzConst=12,
                                           device='cpu')

g_reg = cilPluginToSIRFFactory.getInstance(Diff4th, 
                                           regularisation_parameter=3.5,
                                           edge_parameter=0.02,
                                           iterations=500,
                                           time_marching_parameter=0.001,
                                           device='cpu')
'''
#g_reg = IndicatorBox(lower=0,upper=1)
#%%


n_of_subsets = 14
A = SubsetOperator(am.get_linear_acquisition_model(), n_of_subsets)
A_norms = [PowerMethodNonsquare(Ai, 10, x0=image.copy()) for Ai in A]
#%%
# increase the norms to allow for inaccuracies in their computation
Ls = [1.05 * L for L in A_norms]

f = [KullbackLeibler(op.sirf2sub(noisy_data), op.sirf2sub(background))
     for op in A]

#%%
recon_noreg = spdhg(f, g_noreg, A, A_norms=Ls)

#%%

# set the number of iterations
epochs = 5
niter = epochs * n_of_subsets


# %%
for i in range(niter):
    print(recon_noreg.iter)
    recon_noreg.update()

#%% does currently not work!
recon_reg = spdhg(f, g_reg, A, A_norms=Ls)

# %%
for i in range(niter):
    print(recon_reg.iter)
    recon_reg.update()

#%%  show result
# imshow3(recon_noreg.x, limits=[-0.3,cmax], title='recon noreg')
# imshow3(recon_reg.x, limits=[-0.3,cmax], title='recon reg')
if EXAMPLE in ['thorax', 'brain']:
    sliceno = recon_noreg.x.as_array().shape[0]/2
    fig = plt.figure()
    plt.subplot(1, 3, 1)
    plt.imshow(image.as_array()[sliceno, :, :], cmap='gray')
    plt.clim(-0.3, cmax)
    plt.title('Simulated Image')
    plt.colorbar(shrink=.6)
    plt.axis('off')
    plt.subplot(1, 3, 2)
    plt.imshow(recon_noreg.x.as_array()[sliceno, :, :], cmap='gray')
    plt.clim(-0.3, cmax)
    plt.title('{} iter SPDHG'.format(niter))
    plt.colorbar(shrink=.6)
    plt.axis('off')
    plt.subplot(1, 3, 3)
    plt.imshow(recon_reg.x.as_array()[sliceno, :, :], cmap='gray')
    plt.clim(-0.3, cmax)
    plt.title('{} iter SPDHG + FGP_TV'.format(niter))
    plt.colorbar(shrink=.6)
    plt.axis('off')
    plt.show()
elif EXAMPLE == 'nema':
    # save max for future displays
    cmax = recon_noreg.x.as_array().max()
    sliceno = recon_noreg.x.as_array().shape[0]/2
    sliceno = 71
    fig = plt.figure()
    plt.subplot(1, 2, 1)
    plt.imshow(recon_noreg.x.as_array()[sliceno, :, :], cmap='gray')
    #plt.clim(-0.3,cmax)
    plt.title('{} iter SPDHG'.format(niter))
    plt.colorbar(shrink=.6)
    plt.axis('off')
    plt.subplot(1, 2, 2)
    plt.imshow(recon_reg.x.as_array()[sliceno, :, :], cmap='gray')
    #plt.clim(-0.3,cmax)
    plt.title('{} iter SPDHG + TGV_TV'.format(niter))
    plt.colorbar(shrink=.6)
    plt.axis('off')
    plt.show()
