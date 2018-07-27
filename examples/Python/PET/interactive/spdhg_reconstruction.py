# -*- coding: utf-8 -*-
###
# TBC!!!
# Demonstration of PET reconstruction with CCP PET-MR Software
#
# This demonstration shows how to use OSEM and implement a
# (simplistic) gradient-descent algorithm using SIRF.
#
# This demo is a 'script', i.e. intended to be run step by step in a
# Python IDE such as spyder. It is organised in 'cells'. spyder displays these
# cells nicely and allows you to run each cell on its own.
#
# First version: 8th of September 2016
# Author: Kris Thielemans
#

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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
import numpy
import matplotlib.pyplot as plt
import os
import shutil
#import scipy
#from scipy import optimize
import pSTIR as pet
# plotting settings
plt.ion() # interactive 'on' such that plots appear during loops
# some handy function definitions
def imshow(image, limits, title=''):
    """Usage: imshow(image, [min,max], title)"""
    plt.title(title)
    bitmap=plt.imshow(image)
    if len(limits)==0:
        limits=[image.min(),image.max()]
                
    plt.clim(limits[0], limits[1])
    plt.colorbar(shrink=.6)
    plt.axis('off');
    return bitmap

def make_positive(image_array):
    """truncate any negatives to zero"""
    image_array[image_array<0] = 0;
    return image_array;

def make_cylindrical_FOV(image):
    """truncate to cylindrical FOV"""
    filter = pet.TruncateToCylinderProcessor()
    filter.apply(image)

# go to directory with input files
# adapt this path to your situation (or start everything in the relevant directory)
os.chdir(pet.petmr_data_path('pet'))
# copy files to working folder and change directory to where the output files are
shutil.rmtree('working_folder/thorax_single_slice',True)
shutil.copytree('thorax_single_slice','working_folder/thorax_single_slice')
os.chdir('working_folder/thorax_single_slice')

# We will first create some simulated data from ground-truth images

# Read in images
image = pet.ImageData('emission.hv');
image_array=image.as_array()*.05
image.fill(image_array);
mu_map = pet.ImageData('attenuation.hv');
mu_map_array=mu_map.as_array();
# bitmap display of images
slice=image_array.shape[0]/2;
plt.figure();
#plt.subplot(1,2,1);
imshow(image_array[slice,:,:,], [], 'emission image');
#plt.subplot(1,2,2);
#imshow(mu_map_array[slice,:,:,], [], 'attenuation image');

# save max for future displays
cmax = image_array.max()*.6

# create acquisition model
am = pet.AcquisitionModelUsingRayTracingMatrix()
am.set_num_tangential_LORs(5)
templ = pet.AcquisitionData('template_sinogram.hs')
am.set_up(templ,image); 
# simulate some data using forward projection
acquired_data=am.forward(image)
acquisition_array = acquired_data.as_array()

# Display bitmaps of a middle sinogram
plt.figure()
imshow(acquisition_array[0,:,:,], [], 'Forward projection');

# close all plots
plt.close('all')

#%% create OSMAPOSL reconstructor
# This implements the Ordered Subsets Maximum A-Posteriori One Step Late
# Since we are not using a penalty, or prior in this example, it
# defaults to using MLEM, but we will modify it to OSEM

import pCIL

data = acquired_data.as_array()
background = 0 * data.copy()

f = pCIL.KullbackLeibler(data, background)
g = pCIL.ZeroFun()

init_image = image.asarray().clone()
z = init_image.copy()
y = data.copy()
A = am

niter = 10

L = pCIL.PowerMethodNonsquare(A, 10)
L *= 1.05
tau = 1 / L
sigma = [1 / L]
   
# Selection function
def fun_select(x):
    return [int(numpy.random.choice(len(A), 1, p=1 / len(A)))]
    
recon = pCIL.spdhg(init_image, y, z, f, g, A, tau, sigma, fun_select)

for i in range(niter):
    print(i)
    recon.update()

#%%  create initial image
# we could just use a uniform image but here we will create a disk with a different
# initial value (this will help the display later on)
init_image=image.clone()
init_image.fill(cmax/4)
make_cylindrical_FOV(init_image)
# display
idata = init_image.as_array()
slice=idata.shape[0]/2;
plt.figure()
imshow(idata[slice,:,:],[0,cmax], 'initial image');

#%% reconstruct the image 
reconstructed_image=init_image.clone()
# set up the reconstructor
recon.set_up(reconstructed_image)
# do actual recon
recon.reconstruct(reconstructed_image)

#%% bitmap display of images
reconstructed_array=reconstructed_image.as_array()

plt.figure();
plt.subplot(1,2,1);
imshow(image_array[slice,:,:,], [0,cmax*1.2],'emission image');
plt.subplot(1,2,2);
imshow(reconstructed_array[slice,:,:,], [0,cmax*1.2], 'reconstructed image');


#%% Generate a noisy realisation of the data
noisy_array=numpy.random.poisson(acquisition_array).astype('float64')
print(' Maximum counts in the data: %d' % noisy_array.max())
# stuff into a new AcquisitionData object
noisy_data = acquired_data.clone()
noisy_data.fill(noisy_array);

#%% Display bitmaps of the middle sinogram
plt.figure()
plt.subplot(1,2,1);
imshow(acquisition_array[slice,:,:,], [0,acquisition_array.max()], 'original');
plt.subplot(1,2,2);
imshow(noisy_array[slice,:,:,], [0,acquisition_array.max()], 'noisy');

#%% reconstruct the noisy data
obj_fun.set_acquisition_data(noisy_data)
# We could save the data to file if we wanted to, but currently we don't.
# recon.set_output_filename_prefix('reconstructedImage_noisydata')

noisy_reconstructed_image=init_image.clone()
recon.reconstruct(noisy_reconstructed_image)
#%% bitmap display of images
noisy_reconstructed_array=noisy_reconstructed_image.as_array()

plt.figure();
plt.subplot(1,2,1);
imshow(reconstructed_array[slice,:,:,], [0,cmax*1.2], 'no noise');
plt.subplot(1,2,2);
imshow(noisy_reconstructed_array[slice,:,:,], [0,cmax*1.2], 'with noise');

#%% run same reconstruction but saving images and objective function values every sub-iteration
num_subiters = 64;
all_osem_images = numpy.ndarray(shape=(num_subiters+1,) + idata.shape );
current_image = init_image.clone()
osem_objective_function_values = [ obj_fun.value(current_image) ]
all_osem_images[0,:,:,:] =  current_image.as_array();
for iter in range(1, num_subiters+1):
    recon.update(current_image);
    
    obj_fun_value = obj_fun.value(current_image);
    osem_objective_function_values.append(obj_fun_value);
    all_osem_images[iter,:,:,:] =  current_image.as_array();
  