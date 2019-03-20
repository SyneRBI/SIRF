# -*- coding: utf-8 -*-
###
# Demonstration of basic SPECT capabilities with SIRF:
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
## Copyright 2015 - 2017, 2019 University College London.
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

#%% Initial imports etc
import numpy
from numpy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import sys
import shutil
# plotting settings
plt.ion() # interactive 'on' such that plots appear during loops
#%% Import SIRF
import sirf.STIR
import sirf.Utilities

# Adapt this path to your situation (or start everything in the relevant directory)
os.chdir(sirf.Utilities.examples_data_path('SPECT'))
os.chdir('working_folder/simple')
#%% Read in data
acquired_data=sirf.STIR.AcquisitionData('simulation.hs')
image=acquired_data.create_uniform_image()

#%% create an acquisition model
acq_model_matrix = sirf.STIR.SPECTUBMatrix();
acq_model_matrix.set_keep_all_views_in_cache(True)
acq_model_matrix.set_resolution_model(1,.05, False)
am = sirf.STIR.AcquisitionModelUsingMatrix(acq_model_matrix)
am.set_up(acquired_data,image); 
#%% create objective function
obj_fun = sirf.STIR.make_Poisson_loglikelihood(acquired_data)
obj_fun.set_acquisition_model(am)

#%% create OSMAPOSL reconstructor
# This implements the Ordered Subsets Maximum A-Posteriori One Step Late
# Since we are not using a penalty, or prior in this example, it
# defaults to using MLEM, but we will modify it to OSEM
recon = sirf.STIR.OSMAPOSLReconstructor()
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
sirf.Utilities.show_2D_array('reconstructed image after 5 sub-iterations',reconstructed_array[slice,:,:,]);

#%% do another set of iterations
recon.reconstruct(reconstructed_image)
reconstructed_array=reconstructed_image.as_array()
sirf.Utilities.show_2D_array('reconstructed image after 10 sub-iterations',reconstructed_array[slice,:,:,]);

#%% forward project the image again with this acquisition model and display
acquired_data = am2.forward(image)
acquisition_array = acquired_data.as_array()
sirf.Utilities.show_3D_array(acquisition_array);

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
