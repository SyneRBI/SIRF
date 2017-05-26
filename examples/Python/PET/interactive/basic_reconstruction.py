# -*- coding: utf-8 -*-
"""
Created on Fri May 26 23:06:58 2017

@author: sirfuser
"""


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
bin_efficiencies_array[:,5:10,:] = 0
bin_efficiencies.fill(bin_efficiencies_array)
# now include this in our acquisition model
am.set_bin_efficiency(bin_efficiencies)
#%% forward project the image again with this acquisition model and display
acquired_data = am.forward(image)
acquisition_array = acquired_data.as_array()
plt.figure()
slice_num=acquisition_array.shape[0]/2;
imshow(acquisition_array[slice_num,:,:,], [], 'Forward projection');
