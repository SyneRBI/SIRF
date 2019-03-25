# -*- coding: utf-8 -*-
###
# Demonstration of basic PET capabilities with SIRF: IO and projections
#
# This demonstration shows how to read images and data, display them. It then
# illustrates how to use an AcquisitionModel to forward and backproject.
#
# This demo is a 'script', i.e. intended to be run step by step in a 
# Python IDE such as spyder. It is organised in 'cells'. spyder displays these
# cells nicely and allows you to run each cell on its own.
#
# # We'll use the Python Animation package for one display. This might not display 
# anything depending on your IDE settings (check the 'backend' settings).
# For instance, in spyder, go to Tools->Preferences->iPython->Graphics and
# set your backend to "automatic". You will have to do this BEFORE you start the
# ipython console (or just restart spyder)
#
# Author: Kris Thielemans
# First version: 8th of September 2016
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
#%% Use the 'PET' prefix for all SIRF functions
# This is done here to explicitly differentiate between SIRF pet functions and 
# anything else.
import pSTIR as pet

#%% First define some handy function definitions
# To make subsequent code cleaner, we have a few functions here. You can ignore
# ignore them when you first see this demo.
# They have (minimal) documentation using Python docstrings such that you 
# can do for instance "help(imshow)"
#
# First a function to display an image
def imshow(image, limits, title=''):
    """Display an image with a colourbar, returning the plot handle. 
    
    Arguments:
    image -- a 2D array of numbers
    limits -- colourscale limits as [min,max]. An empty [] uses the full range
    title -- a string for the title of the plot (default "")
    """
    plt.title(title)
    bitmap=plt.imshow(image)
    if len(limits)==0:
        limits=[image.min(),image.max()]
                
    plt.clim(limits[0], limits[1])
    plt.colorbar(shrink=.6)
    plt.axis('off');
    return bitmap

def make_positive(image_array):
    """Truncate any negatives in an ndarray to zero."""
    image_array[image_array<0] = 0;
    return image_array;

def make_cylindrical_FOV(image):
    """Truncate a pet image to a cylindrical FOV."""
    filter = pet.TruncateToCylinderProcessor()
    filter.apply(image)

#%% Go to directory with input files
# Adapt this path to your situation (or start everything in the relevant directory)
os.chdir(pet.examples_data_path('PET'))
#%% Copy files to a working folder and change directory to where these files are.
# We do this to avoid cluttering your SIRF files. This way, you can delete 
# working_folder and start from scratch.
shutil.rmtree('working_folder/brain',True)
shutil.copytree('brain','working_folder/brain')
os.chdir('working_folder/brain')
#%% OK. finally done with initial set-up...

#%% Read in images
# Here we will read some images provided with the demo using the ImageData class.
# These are in Interfile format. Check the main SIRF doc.
image = pet.ImageData('emission.hv');
mu_map = pet.ImageData('attenuation.hv');
#%% What is an ImageData?
# Images are represented by objects with several methods. The most important method 
# is as_array() which we'll use below.
# Let's see what all the methods are.
help(pet.ImageData)

#%% Use as_array to get the underlying array of numbers
image_array=image.as_array();
# This a standard numpy 3D array with its associated methods.
print(image_array.shape)
# Whenever we want to do something with the image-values, we have to do it via this array.
# Let's print a voxel-value.
print(image_array[0,10,20])

#%% Manipulate the image data for illustration
# Multiply the data with a factor
image_array*=0.01;
# Stick this new data into the original image object.
# (This will not modify the file content, only the variable in memory.)
image.fill(image_array);

#%% Display the middle slice of the image (which is really a 3D volume)
# We will use our own imshow function (which was defined above) for brevity.

# Get the middle slice number
slice_num=image_array.shape[0]//2;
# Create a new figure
plt.figure();
# Display the slice
imshow(image_array[slice_num,:,:,], [], 'emission image');

#%% OK. Now we will do some PET projections!
# SIRF uses AcquisitionModel as the object to do forward and back-projections.
# We will create an AcquisitionModel object and then use it to forward project
# our image etc

#%% Create a SIRF acquisition model
# We will use the ray-tracing matrix here as our simple PET model.
# There is more to the accquisition model, but that's for another demo.
am = pet.AcquisitionModelUsingRayTracingMatrix()
# Ask STIR to use 5 LORs per sinogram-element
am.set_num_tangential_LORs(5);

#%% Specify sinogram dimensions
# We need to say what scanner to use, what dimensions etc.
# You do this by using existing PET data as a 'template'. 
# We read a file supplied with the demo as an AcquisitionData object
templ = pet.AcquisitionData('template_sinogram.hs');
# Now set-up our acquisition model with all information that it needs about the data and image.
am.set_up(templ,image); 
#%% The AcquisitionModel is now ready for use

#%% Do a forward projection of our image
# 'forward projection' is the terminology used in PET to simulate the acquisition.
# Input is a SIRF ImageData object (not image_array), output is an AcquisitionData object.
acquired_data=am.forward(image)
#%% Check what methods an AcquisitionData object has
help(acquired_data)
#%% Let's get the Python array
acquisition_array = acquired_data.as_array()
print(acquisition_array.shape)

#%% Display bitmap of the middle sinogram
# AcquisitionData are organised by sinograms, so we need to use the first index
# of the accquisition_array.
plt.figure()
slice_num=acquisition_array.shape[0]//2;
imshow(acquisition_array[slice_num,:,:,], [], 'Forward projection');

#%% Display some different 'views' in a movie
# See note at start of file about your backend if this doesn't work.
bitmaps=[]
fig=plt.figure()
# views are the second index in the data
num_views=acquisition_array.shape[1]
# first construct all the plots
for view in range(0,num_views,4):
    bitmap=plt.imshow(acquisition_array[:,view,:,]);
    plt.clim(0,acquisition_array.max())
    plt.axis('off');
    bitmaps.append([bitmap])
# Display as animation
ani = animation.ArtistAnimation(fig, bitmaps, interval=100, blit=True, repeat_delay=1000);

#%% Let's do a back-projection
# Backprojection uses the transpose of the forward-projection matrix to
# go from AcquisitionData to an ImageData
backprojected = am.backward(acquired_data);
# let's display a slice
plt.figure()
backprojected_array=backprojected.as_array();
imshow(backprojected_array[slice_num,:,:],[], 'backprojection');

#%% close all plots
plt.close('all')

#%% End of this demo!
