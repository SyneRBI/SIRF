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

from pathlib import Path
import numpy as np
import time


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


def main():

	fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
	input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/pDynamicSimulation/'
	output_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'
	
	fpath_xml = input_fpath_prefix + 'Cube128/XCAT_TissueParameters_XML.xml'
	fpath_template_contrast_rawdata = input_fpath_prefix + 'Cube128/CV_nav_cart_128Cube_FLASH_T1.h5'
	
	cartesian_template = False

	if cartesian_template:
		fpath_template_acquisition_rawdata = input_fpath_prefix + 'General/meas_MID29_cart_ref_image_FID78804_ismrmrd.h5'
	else:
		fpath_template_acquisition_rawdata = input_fpath_prefix + 'General/meas_MID30_rad_2d_uniform_FID78805_ismrmrd.h5'
	
	acquisition_ad = pMR.AcquisitionData(fpath_template_acquisition_rawdata)
	if cartesian_template:
		acquisition_ad = pMR.preprocess_acquisition_data(acquisition_ad)
	else:
		acquisition_ad = pMR.set_radial2D_trajectory(acquisition_ad)

	# configure the simulation
	contrast_ad = pMR.AcquisitionData(fpath_template_contrast_rawdata)
	contrast_ad = pMR.preprocess_acquisition_data(contrast_ad)

	labels = pReg.NiftiImageData3D( input_fpath_prefix + "Cube128/label_volume.nii"	)
	mrsim = pDS.MRDynamicSimulation(labels, fpath_xml)

	mrsim.set_contrast_template_data(contrast_ad)
	mrsim.set_acquisition_template_data(acquisition_ad)


	offset_z_mm = -128
	translation = np.array([0, 0, offset_z_mm])
	euler_angles_deg = np.array([15,15,0])

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
	# RESP
	num_sim_resp_states = 2
	resp_motion = pDS.MRMotionDynamic( num_sim_resp_states )
	resp_motion.set_dynamic_signal(t_resp, sig_resp)
	resp_motion.set_cyclicality(False)
	resp_motion.set_groundtruth_folder_prefix(output_fpath_prefix + "output_example_cartesian_3D_simulation_gt_resp")		
	set_motionfields_from_path(resp_motion, input_fpath_prefix + 'Cube128/mvf_resp/')
	mrsim.add_motion_dynamic(resp_motion)

	# CARD
	num_sim_card_states = 2

	card_motion = pDS.MRMotionDynamic(num_sim_card_states)
	card_motion.set_dynamic_signal(t_card, sig_card)
	card_motion.set_cyclicality(True)
	card_motion.set_groundtruth_folder_prefix(output_fpath_prefix + "output_example_cartesian_3D_simulation_gt_card")		
	set_motionfields_from_path(card_motion, input_fpath_prefix + 'Cube128/mvf_card/')
	mrsim.add_motion_dynamic(card_motion)

	#
	fname_sim_output = output_fpath_prefix + "output_example_cartesian_3D_simulation_r{}_c{}.h5".format(num_sim_resp_states,num_sim_card_states)
	simulated_file = Path(fname_sim_output)
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
	recon_nii.write(output_fpath_prefix + "output_example_cartesian_3D_recon_r{}_c{}.nii".format(num_sim_resp_states,num_sim_card_states))
	
	return 1

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
