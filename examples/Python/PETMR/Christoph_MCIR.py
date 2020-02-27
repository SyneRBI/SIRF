import sirf.Gadgetron as pMR
import sirf.STIR as pet
import sirf.Reg as reg
from sirf.Utilities import examples_data_path
from ccpi.utilities.display import plotter2D

pet.AcquisitionData.set_storage_scheme('memory')

#%% Go to directory with input files
# Adapt this path to your situation (or start everything in the relevant directory)
#os.chdir(examples_data_path('PET'))

# import further modules
import os
import numpy as np
# import hdf5storage

import matplotlib.pyplot as plt

import glob

plt.close('all')

# path- and filename of raw data file
#pname = '/media/sf_Data/MRI/Output/PSMR_abstract/'
pname = os.path.abspath(examples_data_path('MR'))
fname = 'grappa2_1rep.h5'


data_path = os.path.join( os.path.abspath('/home/ofn77899') ,
                          'brainweb', 'brainweb_single_slice')

T1_files   = sorted( glob.glob(os.path.join( data_path, 'T1_mf*.nii') ) )
T2_files   = sorted( glob.glob(os.path.join( data_path, 'T2_mf*.nii') ) )
transform_matrices_files  = sorted( glob.glob(os.path.join( data_path, 'fwd_tm*.txt') ) )
for f in T1_files:
    print (f)


# chdir to data dir
# os.chdir(data_path)
# Number of motion states
num_ms = 4

## Create motion states

# Split k-space data into different motion states
# acq = pMR.AcquisitionData(os.path.join(data_path , 'CSM_GRAPPA2_48_FOV180.h5' ))

acq = pMR.AcquisitionData(os.path.join(data_path , 'CSM_FULLY_FOV180.h5' ))

# acq = pMR.AcquisitionData(os.path.join(os.path.abspath('/home/ofn77899') ,
#                           'brainweb' , 'CSM_GRAPPA2_48_FOV180.h5' ) )

# acq = pMR.AcquisitionData(os.path.join(pname , fname ))
acq.sort_by_time()

# Create interleaved sampling
mvec = []
for ind in range(num_ms):
    mvec.append(np.arange(ind, acq.number(), num_ms))

# Go through motion states and create k-space
acq_ms = [0]*num_ms
for ind in range(num_ms):

    acq_ms[ind] = acq.new_acquisition_data(empty=True)

    # Set first two (??) acquisition
    acq_ms[ind].append_acquisition(acq.acquisition(0))
    acq_ms[ind].append_acquisition(acq.acquisition(1))

    # Add motion resolved data
    for jnd in range(len(mvec[ind])):
        if mvec[ind][jnd] < acq.number() - 1 and mvec[ind][jnd] > 1:  # Ensure first and last are not added twice
            cacq = acq.acquisition(mvec[ind][jnd])
            acq_ms[ind].append_acquisition(cacq)

    # Set last acquisition
    acq_ms[ind].append_acquisition(acq.acquisition(acq.number() - 1))

    #acq_ms[ind].write(pname + fname +  '_ms' + str(ind) + '.h5')
    print('MS {} with {} points'.format(ind, acq_ms[ind].number()))


## Calculate reference image and coil maps
prep_data = pMR.preprocess_acquisition_data(acq)

recon = pMR.CartesianGRAPPAReconstructor()
recon.set_input(prep_data)
recon.compute_gfactors(False)
recon.process()
ref_im = recon.get_output()

print ("ref_im shape", ref_im.shape)

fig, ax = plt.subplots(1)
ax.imshow(np.abs(ref_im.as_array()[0,:,:]))

csm = pMR.CoilSensitivityData()
csm.smoothness = 500
csm.calculate(prep_data)

fig, ax = plt.subplots(1,4)
for ind in range(4):
    ax[ind].imshow(np.abs(csm.as_array()[ind,0,:,:]))
plt.show()

t1s = []
for t1f in T1_files:
    t1s.append(
        reg.NiftiImageData(t1f).as_array()
    )
    print (t1s[-1].shape)
plotter2D(t1s)

## Simulate different motion states
ref_im_ms = ref_im.clone()
acq_ms_sim = [0]*num_ms
for ind in range(num_ms):
    AcqMod = pMR.AcquisitionModel(acq_ms[ind], ref_im)
    AcqMod.set_coil_sensitivity_maps(csm)
    ## @EDO: Here you will need your rotation transformation and the BrainWeb image

    if True:
        # The saved brainweb images are real, hence requires casting to complex
        t1 = np.asarray( reg.NiftiImageData(T1_files[ind]).as_array(), dtype=np.complex64)
        # The AcquisitionModel range has shape 1,256,256
        # We plan to put in the brainweb data that are 150,150, so requires this requires padding
        cim = np.pad(t1, ((0,256-150), (0,256-150)), 'constant', constant_values=0)
        # and finally reshape to (1,256,256)
        cim = np.reshape(cim, (1,256,256))
    else:
        cim = ref_im.as_array()
        cim = np.roll(cim, ind*5, axis=1) 
    print ("type", cim.dtype, cim.shape)
    ref_im_ms.fill(cim)

    acq_ms_sim[ind] = AcqMod.forward(ref_im_ms)


## Reconstruct motion states
fig, ax = plt.subplots(1,num_ms+1)

im_sum = 0
AcqModMs = [0]*num_ms
for ind in range(len(acq_ms)):
    prep_data = pMR.preprocess_acquisition_data(acq_ms_sim[ind])

    ## @EDO: With these two lines you define the acquisition model for the individual motion states
    AcqModMs[ind] = pMR.AcquisitionModel(acq_ms[ind], ref_im)
    AcqModMs[ind].set_coil_sensitivity_maps(csm)

    recon_img = AcqModMs[ind].backward(prep_data)
    im = recon_img.as_array()
    im_sum += im
    ax[ind].imshow(np.abs(im[0,:,:]))

ax[num_ms].imshow(np.abs(im_sum[0,:,:]))

plt.show()