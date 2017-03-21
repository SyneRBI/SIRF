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
sys.path.append(os.environ.get('SRC_PATH') + '/SIRF/src/xSTIR/pSTIR')
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
    filter = pet.TruncateToCylinderProcessor()
    filter.apply(image)

#%% go to directory with input files
# adapt this path to your situation (or start everything in the relevant directory)
os.chdir(os.getenv('SRC_PATH')+'/SIRF/src/xSTIR/examples/interactive')
#%% copy files to working folder and change directory to where the output files are
shutil.rmtree('working_folder/single_slice',True)
shutil.copytree('EX_single_slice','working_folder/single_slice')
os.chdir('working_folder/single_slice')
#%% Read in images
image = pet.ImageData('emission.hv');
image_array=image.as_array()*.05
image.fill(image_array);
mu_map = pet.ImageData('CTAC.hv');
mu_map_array=mu_map.as_array();
#%% bitmap display of images
slice=image_array.shape[0]/2;
plt.figure();
#plt.subplot(1,2,1);
imshow(image_array[slice,:,:,], [], 'emission image');
#plt.subplot(1,2,2);
#imshow(mu_map_array[slice,:,:,], [], 'attenuation image');

#%% save max for future displays
cmax = image_array.max()*.6

#%% create matrix to be used by the acquisition model
matrix = pet.RayTracingMatrix()
matrix.set_num_tangential_LORs(5)

#%% create acquisition model
am = pet.AcquisitionModelUsingMatrix()
am.set_matrix(matrix)
templ = pet.AcquisitionData('template_sinogram.hs')
am.set_up(templ,image); 
#%% do forward projection
acquired_data=am.forward(image)
acquisition_array = acquired_data.as_array()

#%% Display bitmaps of a middle sinogram
plt.figure()
imshow(acquisition_array[0,:,:,], [], 'Forward projection');

#%% close all plots
plt.close('all')

#%% create objective function
obj_fun = pet.PoissonLogLh_LinModMean_AcqMod()
obj_fun.set_acquisition_model(am)
obj_fun.set_acquisition_data(acquired_data)
#obj_fun.set_prior(prior)

#%% create OSMAPOSL reconstructor
# This defaults to using EM, but we will modify it to OSEM
recon = pet.OSMAPOSLReconstructor()
recon.set_objective_function(obj_fun)
recon.set_num_subsets(2)
num_iters=10;
recon.set_num_subiterations(num_iters)
recon.set_output_filename_prefix('reconstructedImage')
#%%  create initial image
# just create a disk
init_image=image.get_empty_copy()
init_image.fill(cmax/4)
make_cylindrical_FOV(init_image)
# display
idata = init_image.as_array()
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

slice=image_array.shape[0]/2;
plt.figure();
plt.subplot(1,2,1);
imshow(image_array[slice,:,:,], [0,cmax*1.2],'emission image');
plt.subplot(1,2,2);
imshow(reconstructed_array[slice,:,:,], [0,cmax*1.2], 'reconstructed image');


#%% Add some noise to the data
noisy_array=numpy.random.poisson(acquisition_array).astype('float64')
print(' Maximum counts in the data: %d' % noisy_array.max())
# stuff into a new AcquisitionData object
noisy_data = acquired_data.get_empty_copy()
noisy_data.fill(noisy_array);

#%% Display bitmaps of the middle sinogram
plt.figure()
plt.subplot(1,2,1);
imshow(acquisition_array[slice,:,:,], [0,acquisition_array.max()], 'original');
plt.subplot(1,2,2);
imshow(noisy_array[slice,:,:,], [0,acquisition_array.max()], 'noisy');

#%% reconstruct the noisy data
obj_fun.set_acquisition_data(noisy_data)
recon.set_output_filename_prefix('reconstructedImage_noisydata')

noisy_reconstructed_image=init_image.clone()
recon.reconstruct(noisy_reconstructed_image)
#%% bitmap display of images

noisy_reconstructed_array=noisy_reconstructed_image.as_array()

plt.figure();
plt.subplot(1,2,1);
imshow(reconstructed_array[slice,:,:,], [0,cmax*1.2], 'no noise');
plt.subplot(1,2,2);
imshow(noisy_reconstructed_array[slice,:,:,], [0,cmax*1.2], 'with noise');

#%% run same reconstruction but saving images and objective function values every iteration
num_iters = 64;
all_em_images = numpy.ndarray(shape=(num_iters+1,) + idata.shape );
current_image = init_image.clone()
em_objective_function_values = [ obj_fun.value(current_image) ]
all_em_images[0,:,:,:] =  current_image.as_array();
for iter in range(1, num_iters+1):
    recon.update(current_image);
    
    obj_fun_value = obj_fun.value(current_image);
    em_objective_function_values.append(obj_fun_value);
    all_em_images[iter,:,:,:] =  current_image.as_array();
  
#%% define a function for plotting images and the updates
def plot_progress(all_images, title, iterations = []):
    if len(iterations)==0:
        num_iters = all_images[0].shape[0]-1;
        iterations = range(1, num_iters+1);
    num_rows = len(all_images);
    plt.close('all');
    for iter in iterations:
        plt.figure(iter)
        for r in range(num_rows):
            plt.subplot(num_rows,2,2*r+1)
            imshow(all_images[r][iter,slice,:,:], [0,cmax], '%s at %d' % (title[r],  iter))
            plt.subplot(num_rows,2,2*r+2)
            imshow(all_images[r][iter,slice,:,:]-all_images[r][iter-1,slice,:,:],[-cmax*.1,cmax*.1], 'update')
            plt.pause(.05)
        
#%% now call this function to see how we went along
iterations = (1,2,4,8,16,32,64);
plot_progress([all_em_images], ['OSEM'],iterations)

#%% plot objective function values
plt.figure()
#plt.plot(iterations, [ em_objective_function_values[i] for i in iterations])
plt.plot(em_objective_function_values)
plt.title('Objective function values')
plt.xlabel('iterations')

#%% ROI
ROI_lesion = all_em_images[:,(slice,), 65:70, 40:45];
ROI_lung = all_em_images[:,(slice,), 75:80, 45:50];

ROI_mean_lesion = ROI_lesion.mean(axis=(1,2,3))
ROI_std_lesion = ROI_lesion.std(axis=(1,2,3))

ROI_mean_lung = ROI_lung.mean(axis=(1,2,3))
ROI_std_lung = ROI_lung.std(axis=(1,2,3))

plt.figure()
plt.hold('on')
plt.subplot(1,2,1)
plt.plot(ROI_mean_lesion,'k')
plt.plot(ROI_mean_lung,'r')
plt.title('ROI mean')
plt.xlabel('iterations')
plt.subplot(1,2,2)
plt.plot(ROI_std_lesion, 'k')
plt.plot(ROI_std_lung, 'r')
plt.title('ROI standard deviation')
plt.xlabel('iterations');

#%%
plt.close('all')









#%% Perform gradient descent for 8 iterations
# Gradient descent (GD) works by updating the image in the direction of the gradient
#    new_image = current_image + step_size * gradient

num_iters = 32
# relative step-size
tau = .3

# set initial image and store it and the value of the objective function for plotting
current_image = init_image.clone()
GD_objective_function_values = [ obj_fun.value(current_image) ]
idata = current_image.as_array()
all_images = numpy.ndarray(shape=(num_iters+1,) + idata.shape );
all_images[0,:,:,:] =  idata;

#%% perform GD iterations
for iter in range(1, num_iters+1):  
    # obtain gradient
    grad = obj_fun.gradient(current_image, 0)
    grad_array = grad.as_array()
    
    # compute step-size as relative to current image-norm
    step_size = tau * norm(idata) / norm(grad_array)
    
    # perform gradient descent step and truncate to positive values
    idata = make_positive(idata + step_size*grad_array)
    current_image.fill(idata)

    # compute objective function value for plotting, and write some diagnostics
    obj_fun_value = obj_fun.value(current_image)
    GD_objective_function_values.append(obj_fun_value)
    all_images[iter,:,:,:] = idata; 
#%% Plot objective function values
plt.figure()
plt.hold('on')
plt.title('Objective function value vs iterations')
plt.plot(GD_objective_function_values,'b');
plt.plot(em_objective_function_values,'r');
plt.legend(('gradient descent', 'OSEM'),loc='lower right');

#%% compare GD and OSEM images
plot_progress([all_images, all_em_images], ['GD' ,'OSEM'],[2,4,8,16,32])
