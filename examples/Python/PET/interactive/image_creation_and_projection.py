# -*- coding: utf-8 -*-
###
# Demonstration of basic PET capabilities with SIRF: 
# creating images using shapes and project them
#
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
#%% create a shape
shape = pet.EllipticCylinder()
# define its size (in mm)
shape.set_length(50)
shape.set_radii((30, 40))
# centre of shape in (x,y,z) coordinates where (0,0,0) is centre of first plane
shape.set_origin((60, -30, 20))

#%% add the shape to the image
# first set the image values to 0
image.fill(0)
image.add_shape(shape, scale = 1)

#%% add same shape at different location and with different intensity
shape.set_origin((-60, -30, 40))
image.add_shape(shape, scale = 0.75)

#%% show the phantom image as a sequence of transverse images
show_3D_array(image.as_array())

#%% forward project this image and display all sinograms
acquired_data = am.forward(image)
acquisition_array = acquired_data.as_array();
show_3D_array(acquisition_array);
#%% Show every 8th view 
# Doing this here with a complicated one-liner...
show_3D_array(acquisition_array[:,range(0,acquisition_array.shape[1],8),:].transpose(1,0,2))
# You could now of course try the animation of the previous demo...
