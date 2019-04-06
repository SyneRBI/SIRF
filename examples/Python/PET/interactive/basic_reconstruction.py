# -*- coding: utf-8 -*-
###
# Demonstration of basic PET capabilities with SIRF: 
# basic OSEM reconstruction, projection with a (slightly) more sophisticated acquisition model
##
# This demo is a 'script', i.e. intended to be run step by step in a 
# Python IDE such as spyder. It is organised in 'cells'. spyder displays these
# cells nicely and allows you to run each cell on its own.
#
# WARNING: This script assumes you have run the display_and_projection.py demo first!
#
# Author: Kris Thielemans
# Author: Evgueni Ovtchinnikov
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

#%% just check if you ran the previous demo
if 'image' in globals():
    print('Ok, we can proceed')
else:
    print('This script assumes you have run the display_and_projection.py demo first!')

#%% Import some extra functions
from pUtilities import show_2D_array, show_3D_array

#%% Do a forward projection of our image
am = pet.AcquisitionModelUsingRayTracingMatrix()
am.set_up(templ,image); 
acquired_data=am.forward(image)
#%% create objective function
obj_fun = pet.make_Poisson_loglikelihood(acquired_data)
obj_fun.set_acquisition_model(am)

#%% create OSMAPOSL reconstructor
# This implements the Ordered Subsets Maximum A-Posteriori One Step Late
# Since we are not using a penalty, or prior in this example, it
# defaults to using MLEM, but we will modify it to OSEM
recon = pet.OSMAPOSLReconstructor()
recon.set_objective_function(obj_fun)
recon.set_num_subsets(4)
recon.set_num_subiterations(5)

#%% reconstruct the image 
# First create a new image to use for the reconstruction
# We will just use the original as a 'template' to have the same voxel sizes etc
reconstructed_image=image.clone()
# Set its values to 1 to create a uniform image
reconstructed_image.fill(1)
# set up the reconstructor
recon.set_up(reconstructed_image)
# do actual recon
recon.reconstruct(reconstructed_image)

#%% display of image
reconstructed_array=reconstructed_image.as_array()
slice=reconstructed_array.shape[0]/3;
show_2D_array('reconstructed image after 5 sub-iterations',reconstructed_array[slice,:,:,]);

#%% do a another set of iterations
recon.reconstruct(reconstructed_image)
reconstructed_array=reconstructed_image.as_array()
show_2D_array('reconstructed image after 10 sub-iterations',reconstructed_array[slice,:,:,]);

#%% We now add a multiplicative term to the acquisition model
# In PET, detector-pairs have different efficiencies. We want to include
# this in our 'forward' model such that the reconstruction can
# take this into account.
#
# The way to do this in SIRF is to include 'bin efficiencies' in the model,
# i.e. one multiplicative factor for each bin in the data.
#
# You would normally derive these efficiencies from a "normalisation" scan.
# Here we will simply set the efficiencies for some 'views' to zero.
# This is actually physically impossible for PET (although ok for SPECT),
# but this is only a demo!

# first create a copy of the data such that we have an object of the appropriate size
bin_efficiencies = acquired_data.clone()
# set all values to 1
bin_efficiencies.fill(1.)
# set a portion of bin efficiencies to zero;
bin_efficiencies_array = bin_efficiencies.as_array()
bin_efficiencies_array[0,:,5:20,:] = 0
bin_efficiencies.fill(bin_efficiencies_array)
#%% Create a new acquisition model
am2 = pet.AcquisitionModelUsingRayTracingMatrix()
am2.set_num_tangential_LORs(5);
am2.set_up(templ,image); 
# now include the bin efficiencies in our acquisition model
asm = pet.AcquisitionSensitivityModel(bin_efficiencies)
am2.set_acquisition_sensitivity(asm)
am2.set_up(templ,image);
#%% forward project the image again with this acquisition model and display
acquired_data = am2.forward(image)
acquisition_array = acquired_data.as_array()
show_3D_array(acquisition_array);

#%% Let us reconstruct this data with the original acquisition model (without bin efficiencies)
obj_fun.set_acquisition_data(acquired_data)
obj_fun.set_acquisition_model(am)
reconstructed_image.fill(1)
recon.set_up(reconstructed_image)
recon.set_num_subiterations(10)
recon.reconstruct(reconstructed_image)
#%% display
# we fix the max for the colour scale related to the true max
cmax = image.as_array().max()*1.2;
reconstructed_array=reconstructed_image.as_array()
plt.figure()
imshow(reconstructed_array[slice,:,:,], [0,cmax],'reconstructed image with original acquisition model');
#%% Now we use the correct acquisition model
obj_fun.set_acquisition_model(am2)
reconstructed_image.fill(1)
recon.set_up(reconstructed_image)
recon.reconstruct(reconstructed_image)
#%% display
reconstructed_array=reconstructed_image.as_array()
plt.figure()
imshow(reconstructed_array[slice,:,:,], [0,cmax],'reconstructed image with correct acquisition model');
