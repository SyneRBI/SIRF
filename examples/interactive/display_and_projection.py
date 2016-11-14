# -*- coding: utf-8 -*-
###
# Demonstration of PET reconstruction with CCP PET-MR Software
# 8th of September 2016
#
# Author: Kris Thielemans
#
#%% Initial imports etc
import numpy
from numpy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import sys
import shutil
#import scipy
#from scipy import optimize
sys.path.append(os.environ.get('SRC_PATH') + '/xSTIR/pSTIR')
import pStir as pet
# plotting settings
plt.ion() # interactive 'on' such that plots appear during loops
#%% some handy function definitions
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
    filter = pet.CylindricFilter()
    filter.apply(image)

#%% go to directory with input files
# adapt this path to your situation (or start everything in the relevant directory)
os.chdir(os.getenv('SRC_PATH')+'/xSTIR/examples/interactive')
#%% copy files to working folder and change directory to where the output files are
shutil.rmtree('working_folder/brain',True)
shutil.copytree('EX_brain','working_folder/brain')
os.chdir('working_folder/brain')
#%% Read in images
image = pet.Image('emission.hv');
image_array=image.as_array()*.01;
image.fill(image_array);
mu_map = pet.Image('attenuation.hv');
mu_map_array=mu_map.as_array();
#%% bitmap display of images
slice=image_array.shape[0]/2;
plt.figure();
#plt.subplot(1,2,1);
imshow(image_array[slice,:,:,], [], 'emission image');
#plt.subplot(1,2,2);
#imshow(mu_map_array[slice,:,:,], [], 'attenuation image');

#%% save max for future displays
cmax = image_array.max()

#%% create matrix to be used by the acquisition model
matrix = pet.RayTracingMatrix()
matrix.set_num_tangential_LORs(5)

#%% create acquisition model
am = pet.AcquisitionModelUsingMatrix()
am.set_matrix(matrix)
templ = pet.AcquisitionData('template_sinogram.hs');
am.set_up(templ,image); 
#%% do forward projection
acquired=am.forward(image)
acquisition_array = acquired.as_array()

#%% Display bitmaps of a middle sinogram
plt.figure()
imshow(acquisition_array[slice,:,:,], [], 'Forward projection');

#%% Display some different views in a movie
bitmaps=[]
fig=plt.figure()
for view in range(0,64,4):
    bitmap=plt.imshow(acquisition_array[:,view,:,]);
    plt.clim(0,acquisition_array.max())
    plt.axis('off');
    bitmaps.append([bitmap])

ani = animation.ArtistAnimation(fig, bitmaps, interval=100, blit=True, repeat_delay=1000);

#%% Let's do a back-projection
backprojected_array = am.backward(acquired).as_array()
plt.figure()
imshow(backprojected_array[slice,:,:],[], 'backprojection');

#%% close all plots
plt.close('all')

