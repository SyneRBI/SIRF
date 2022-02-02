'''
bla bla

Usage:
cartesian_3D_simulation.py [--help | options]

Options:
--non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF).
## Copyright 2015 - 2020 Rutherford Appleton Laboratory STFC.
## Copyright 2015 - 2017 University College London.
## Copyright 2015 - 2017 Physikalisch-Technische Bundesanstalt.
##
## This is software developed for the Collaborative Computational
## Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
## (http://www.ccpsynerbi.ac.uk/).
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

__version__ = '0.1.0'
from docopt import docopt

args = docopt(__doc__, version=__version__)

from pUtilities import *
import sirf.Reg as pReg
import sirf.DynamicSimulation as pDS
import sirf.Gadgetron as pMR

# import engine module

# process command-line options
show_plot = not args['--non-interactive']

import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import time
import random

def create_dummy_mrf_signal(num_tissue_types,num_time_points):
    # 
    labels = np.array([i for i in range(num_tissue_types)])
    time_curve = np.sin(np.array([t for t in range(num_time_points)]) * np.pi / num_time_points) \
                +1j * np.cos(np.array([t for t in range(num_time_points)]) * 2 * np.pi / num_time_points)

    mrf_signal = labels[:,np.newaxis] * time_curve[np.newaxis,:] / num_tissue_types / np.sqrt(2)
    
    return pDS.ExternalMRSignal(labels, mrf_signal)

def read_motionfields(fpath_prefix):
	p = sorted( Path(fpath_prefix).glob('*.nii') )
	files = [x for x in p if x.is_file()]
	
	temp = []
	for f in files:
		print("Reading from {} ... ".format(f))
		img = pReg.NiftiImageData3DDisplacement(str(f))
		temp.append(img)

	data = np.array(temp, dtype=object)
	return data

def set_motionfields_from_path(modyn, fpath_prefix):

	assert_validity(modyn, pDS.MRMotionDynamic)
	mvfs = read_motionfields(fpath_prefix)

	for m in mvfs:
		modyn.add_displacement_field(m)


def get_normed_surrogate_signal(t0_s, tmax_s, Nt, f_Hz):

	t_s = np.linspace(t0_s, tmax_s, Nt)
	sig = 0.5 * (1 + np.sin( 2*np.pi*f_Hz*t_s))
	return t_s, sig


def static_MR_fingerprinting():

    print(" --- Running static MRF simulation.\n")

    fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
    prefix_fpath_input = fpath_testdata_prefix + 'Input/xDynamicSimulation/pDynamicSimulation/'
    
    fpath_output = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'
    fname_output = fpath_output + 'simulated_static_MR_fingerprinting.h5'

    fname_xml = prefix_fpath_input + 'Cube128/XCAT_TissueParameters_XML.xml'
    fname_template_contrast_rawdata = prefix_fpath_input + 'Cube128/CV_nav_cart_128Cube_FLASH_T1.h5'
    fname_template_acquisition_rawdata = prefix_fpath_input + 'General/meas_MID27_CV_11s_TI2153_a6_2x2x8_TR45_FID33312.h5'

    print(" --- Loading the template raw data.\n")

    contrast_template = pMR.AcquisitionData(fname_template_contrast_rawdata)
    contrast_template = pMR.preprocess_acquisition_data(contrast_template)

    acquisition_template = pMR.AcquisitionData(fname_template_acquisition_rawdata)
    acquisition_template = pMR.set_radial2D_trajectory(acquisition_template)
    

    # First compute the coil profiles from the rawdata itself.
    csm = pMR.CoilSensitivityData()
    csm.calculate(acquisition_template)


    print(" --- Extracting subset for template acquisition data.\n")
    num_acquisitions = 3
    list_acquisitions_to_keep = np.array(random.sample(range(0,acquisition_template.number()), num_acquisitions))
    print("We will keep {} acquisitions.".format(list_acquisitions_to_keep.size))

    simulated_traj = pMR.get_data_trajectory(acquisition_template)
    plt.figure()
    plt.scatter(simulated_traj[:,0], simulated_traj[:,1] )
    plt.title('prior')
    plt.show()

    reduced_template = acquisition_template.get_subset(list_acquisitions_to_keep) 
    reduced_template = pMR.set_goldenangle2D_trajectory(reduced_template)
    simulated_traj = pMR.get_data_trajectory(reduced_template)
    plt.figure('posterior')
    plt.scatter(simulated_traj[:,0], simulated_traj[:,1] )
    plt.show()

    print(" --- Our Acquisition template contains {} readouts.\n".format(reduced_template.number()))
    ##
    labels = pReg.NiftiImageData3D( prefix_fpath_input + "Cube128/label_volume.nii"	)
    mrsim = pDS.MRDynamicSimulation(labels, fname_xml)

    mrsim.set_csm(csm)
    mrsim.set_contrast_template_data(contrast_template)
    mrsim.set_acquisition_template_data(acquisition_template)

    offset_z_mm = -128
    offset_centre_mm = 0
    translation = np.array([offset_centre_mm, offset_centre_mm, offset_z_mm])
    euler_angles_deg = np.array([0,0,0])

    offset_trafo = pReg.AffineTransformation(translation, euler_angles_deg)
    mrsim.set_offset_trafo(offset_trafo)

    # set which tissue defines SNR
    SNR = 10
    SNR_label = 13

    mrsim.set_snr(SNR)
    mrsim.set_snr_label(SNR_label)

    fname_simulation_output = "simulation_static_MRF"
    simulated_file = Path(fpath_output, fname_simulation_output).with_suffix('.h5')
    if not simulated_file.is_file():
        print(" --- Simulation of MRF data \n")
        tstart = time.time()
        mrsim.simulate_data()
        print("--- Required {} minutes for the simulation.".format( (time.time()-tstart)/60))
        mrsim.write_simulation_results(str(simulated_file))
    else:
        print("Skipping simulation since output file already exists.")

    simulated_data = pMR.AcquisitionData(str(simulated_file))
    csm.calculate(simulated_data)

    AM = pMR.AcquisitionModel()
    AM.set_coil_sensitivity_maps(csm)
    AM.set_up(simulated_data, csm)

    recon_img = AM.inverse(simulated_data)
    recon_nii = pReg.NiftiImageData3D(recon_img)
    recon_nii = recon_nii.abs()

    fname_output = fpath_output + "recon_" + fname_simulation_output + ".nii"
    recon_nii.write(fname_output)

    return 1

def motion_MR_fingerprinting():

    fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
    input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/pDynamicSimulation/'
    output_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'

    fpath_xml = input_fpath_prefix + 'Slab128/XCAT_TissueParameters_XML.xml'
    fpath_template_contrast_rawdata = input_fpath_prefix + 'Slab128/CV_nav_cart_128Slab_FLASH_T1.h5'

    trajectory_type = 'radial2D'

    if trajectory_type == 'cartesian':
        fpath_template_acquisition_rawdata = input_fpath_prefix + 'General/meas_MID29_cart_ref_image_FID78804_ismrmrd.h5'
    else:
        fpath_template_acquisition_rawdata = input_fpath_prefix + 'General/meas_MID30_rad_2d_uniform_FID78805_ismrmrd.h5'

    acquisition_ad = pMR.AcquisitionData(fpath_template_acquisition_rawdata)

    if trajectory_type == 'cartesian':
        acquisition_ad = pMR.preprocess_acquisition_data(acquisition_ad)
    elif trajectory_type == 'radial2D':
        acquisition_ad = pMR.set_radial2D_trajectory(acquisition_ad)
    elif trajectory_type == 'goldenangle2D':
        acquisition_ad = pMR.set_goldenangle2D_trajectory(acquisition_ad)
    else:
        raise ValueError("The trajectory you gave is {}. The only options are cartesian, radial2D or goldenangle2D".format(trajectory_type))

    # configure the simulation
    contrast_ad = pMR.AcquisitionData(fpath_template_contrast_rawdata)
    contrast_ad = pMR.preprocess_acquisition_data(contrast_ad)

    labels = pReg.NiftiImageData3D( input_fpath_prefix + "Slab128/label_volume.nii"	)
    mrsim = pDS.MRDynamicSimulation(labels, fpath_xml)

    mrsim.set_contrast_template_data(contrast_ad)
    mrsim.set_acquisition_template_data(acquisition_ad)

    offset_z_mm = 0
    translation = np.array([64, 64, offset_z_mm])
    euler_angles_deg = np.array([0,3,3])

    offset_trafo = pReg.AffineTransformation(translation, euler_angles_deg)
    mrsim.set_offset_trafo(offset_trafo)

    # take CSM from the rawdata itself
    # could be replaced if independent way of computing CSM is available
    csm = pMR.CoilSensitivityData()
    csm.calculate(acquisition_ad)
    mrsim.set_csm(csm)

    # set which tissue defines SNR
    SNR = 10
    SNR_label = 13

    mrsim.set_snr(SNR)
    mrsim.set_snr_label(SNR_label)

    # configure the surrogates
    Nt = 10000
    t0_s = 0
    tmax_s = 60* 5 

    f_Hz_card = 1
    f_Hz_resp = 0.2

    t_resp, sig_resp = get_normed_surrogate_signal(t0_s, tmax_s, Nt, f_Hz_resp)
    t_card, sig_card = get_normed_surrogate_signal(t0_s, tmax_s, Nt, f_Hz_card)

    # configure the motion
    num_motion_states = 2
    # RESP
    num_sim_resp_states = num_motion_states
    resp_motion = pDS.MRMotionDynamic( num_sim_resp_states )
    resp_motion.set_dynamic_signal(t_resp, sig_resp)
    resp_motion.set_cyclicality(False)
    resp_motion.set_groundtruth_folder_prefix(output_fpath_prefix + "output_simulation_2D_motiondata_r_{}_gt_resp".format(num_sim_resp_states))		
    set_motionfields_from_path(resp_motion, input_fpath_prefix + 'Slab128/mvf_resp/')
    mrsim.add_motion_dynamic(resp_motion)

    # external signal
    num_max_labels = 100
    external_signal = create_dummy_mrf_signal(num_max_labels, acquisition_ad.number())

    external_contrast = pDS.ExternalMRContrastDynamic() 
    external_contrast.add_external_signal(external_signal)

    mrsim.add_external_contrast_dynamic(external_contrast)


    # CARD
    num_sim_card_states = num_motion_states

    card_motion = pDS.MRMotionDynamic(num_sim_card_states)
    card_motion.set_dynamic_signal(t_card, sig_card)
    card_motion.set_cyclicality(True)
    card_motion.set_groundtruth_folder_prefix(output_fpath_prefix + "output_simulation_2D_externalcontrast_c_{}_gt_card".format(num_sim_card_states))		
    set_motionfields_from_path(card_motion, input_fpath_prefix + 'Slab128/mvf_card/')
    mrsim.add_motion_dynamic(card_motion)

    #
    fname_simulation_output = "output_simulation_2D_externalcontrast_traj_{}_r{}_c{}".format(trajectory_type, num_sim_resp_states,num_sim_card_states)

    fname_output = output_fpath_prefix + fname_simulation_output + ".h5"
    simulated_file = Path(fname_output)
    if not simulated_file.is_file():

        tstart = time.time()
        mrsim.simulate_data()
        print("--- Required {} minutes for the simulation.".format( (time.time()-tstart)/60))
        mrsim.write_simulation_results(str(simulated_file))
    else:
        print("Skipping simulation since output file already exists.")

    mrsim.save_motion_ground_truth()

    simulated_data = pMR.AcquisitionData(str(simulated_file))

    csm.calculate(simulated_data)

    AM = pMR.AcquisitionModel()
    AM.set_coil_sensitivity_maps(csm)
    AM.set_up(simulated_data, csm)

    recon_img = AM.inverse(simulated_data)
    recon_nii = pReg.NiftiImageData3D(recon_img)
    recon_nii = recon_nii.abs()

    fname_output = output_fpath_prefix + "recon_" + fname_simulation_output + ".nii"
    recon_nii.write(fname_output)

    return 1

def main():
    static_MR_fingerprinting()
    # motion_MR_fingerprinting()
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
