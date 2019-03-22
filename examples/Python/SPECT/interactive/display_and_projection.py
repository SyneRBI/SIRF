# -*- coding: utf-8 -*-
###
# Demonstration of basic SPECT capabilities with SIRF: IO and projections
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
#%% import the SIRF module into Python
import sirf.Utilities
import sirf.STIR
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

def create_sample_image(image):
    '''
    fill the image with some simple geometric shapes
    '''
    image.fill(0)
    # create a shape
    shape = sirf.STIR.EllipticCylinder()
    shape.set_length(400)
    shape.set_radii((40, 100))
    shape.set_origin((10, 60, 0))

    # add the shape to the image
    image.add_shape(shape, scale = 1)

    # add another shape
    shape.set_radii((30, 30))
    shape.set_origin((10, -30, 60))
    image.add_shape(shape, scale = 1.5)

    # add another shape
    shape.set_origin((10, -30, -60))
    image.add_shape(shape, scale = 0.75)

    # another
    shape = sirf.STIR.EllipticCylinder()
    shape.set_length(40)
    shape.set_radii((20, 20))
    shape.set_origin((30, 30, 0))
    image.add_shape(shape, scale = 1)

def create_attenuation_image(image):
    '''
    fill the attenuation image with some simple geometric shapes
    '''
    image.fill(0)
    # create a shape
    shape = sirf.STIR.EllipticCylinder()
    shape.set_length(400)
    shape.set_radii((120, 120))
    shape.set_origin((0, 0, 0))

    # add the shape to the image
    image.add_shape(shape, scale = 0.15)

    # add another shape
    shape.set_radii((30, 30))
    shape.set_origin((60, -30, 10))
    image.add_shape(shape, scale = -.1)

#%% Go to directory with input files
# Adapt this path to your situation (or start everything in the relevant directory)
os.chdir(sirf.Utilities.examples_data_path('SPECT'))
#%% Copy files to a working folder and change directory to where these files are.
# We do this to avoid cluttering your SIRF files. This way, you can delete
# working_folder and start from scratch.
dest_dir='working_folder/simple';
shutil.rmtree(dest_dir,True)
os.makedirs(dest_dir)
import glob
for file in glob.glob(r'*s'): # copy all *s files
    shutil.copy(file, dest_dir)
os.chdir('working_folder/simple')

#%% OK. finally done with initial set-up...

#%% Specify sinogram dimensions
# We need to say what scanner to use, what dimensions etc.
# You do this by using existing spect data as a 'template'.
# We read a file supplied with the demo as an AcquisitionData object
acq_template = sirf.STIR.AcquisitionData('template_sinogram.hs');
acq_dimensions = acq_template.dimensions()
#%% Create an emission and attenuation image with suitable sizes
image = sirf.STIR.ImageData();
image_size = (acq_dimensions[0],acq_dimensions[2], acq_dimensions[2]);
voxel_size = (3.32, 3.32, 3.32); # TODO: get these from acq_template
image.initialise(image_size, voxel_size)
create_sample_image(image)
image.write("emission.hv")

atten_image = image.get_uniform_copy(0)
create_attenuation_image(atten_image)
atten_image.write("attenuation.hv")

#%% What is an ImageData?
# Images are represented by objects with several methods. The most important method
# is as_array() which we'll use below.
# Let's see what all the methods are.
help(sirf.STIR.ImageData)

#%% Use as_array to get the underlying array of numbers
image_array=image.as_array();
# This a standard numpy 3D array with its associated methods.
print(image_array.shape)
# Whenever we want to do something with the image-values, we have to do it via this array.
# Let's print a voxel-value.
print(image_array[0,80,60])

#%% Manipulate the image data for illustration
# Multiply the data with a factor
image_array*=0.01;
# Stick this new data into the original image object.
# (This will not modify the file content, only the variable in memory.)
image.fill(image_array);

#%% Display the middle slice of the images (which is really a 3D volume)
# We will use our own imshow function (which was defined above) for brevity.

# Get the middle slice number
slice_num=image_array.shape[0]//2;
# Create a new figure
plt.figure();
# Display the slice
plt.subplot(1,2,1)
imshow(image_array[slice_num,:,:,], [], 'emission image');
plt.subplot(1,2,2)
imshow(atten_image.as_array()[slice_num,:,:,], [], 'attenuation image');
#%% Or show all slices
sirf.Utilities.show_3D_array(image_array,suptitle='emission image')
#%% OK. Now we will do some SPECT projections!
# SIRF uses AcquisitionModel as the object to do forward and back-projections.
# We will create an AcquisitionModel object and then use it to forward project
# our image etc

#%% Create a SIRF acquisition model
# We use a matrix developed by people from the University of Barvelona as our system model.
acq_model_matrix = sirf.STIR.SPECTUBMatrix();

#%% Set some options
# you can find out more about these by typing asking for help
help(acq_model_matrix.set_resolution_model)
acq_model_matrix.set_resolution_model(1,.05, False)

#%% Add the attenuation to the model
acq_model_matrix.set_attenuation_image(atten_image)
#%% create an acquisition model with this matrix
am = sirf.STIR.AcquisitionModelUsingMatrix(acq_model_matrix)
#%% Now set-up our acquisition model with all information that it needs about the data and image.
am.set_up(acq_template,image);
#%% The AcquisitionModel is now ready for use

#%% Do a forward projection of our image
# 'forward projection' is the terminology used in spect to simulate the acquisition.
# Input is a SIRF ImageData object (not image_array), output is an AcquisitionData object.
acquired_data=am.forward(image)
#%% Check what methods an AcquisitionData object has
help(acquired_data)
#%% Write it to file for future use
acquired_data.write('simulation.hs')
#%% Let's get the Python array
acquisition_array = acquired_data.as_array()
print(acquisition_array.shape)

#%% Display bitmap of the middle sinogram
# AcquisitionData are organised by sinograms, so we need to use the first index
# of the accquisition_array.
plt.figure()
slice_num=acquisition_array.shape[0]//2;
imshow(acquisition_array[slice_num,:,:,], [], 'Forward projection');
#%% Show all sinograms
sirf.Utilities.show_3D_array(acquisition_array, suptitle='Forward projection');

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
