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
from numpy.core.numeric import identity

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


def offset_shift_from_nifti(float_img):


    # rot_angle_rad = np.pi/4
    # ca = np.cos(rot_angle_rad)
    # sa = np.sin(rot_angle_rad)

    # tm_affine = np.array([[ca,-sa,0,0],
    #                       [sa,ca,0,0],
    #                       [0,0   ,1,0],
    #                       [0,0   ,0,1]])
    
    # return pReg.AffineTransformation(tm_affine)

    shift_array = -1 * float_img.get_voxel_sizes()[1:4] * np.array(float_img.shape)/2
    print("--- Resolution is {} mm".format(float_img.get_voxel_sizes()))
    print("--- We will shift by {} mm".format(shift_array))

    shift_x = pReg.NiftiImageData3D(float_img)
    shift_y = pReg.NiftiImageData3D(float_img)
    shift_z = pReg.NiftiImageData3D(float_img)

    shift_x.fill(0)
    shift_y.fill(0)
    shift_z.fill(shift_array[2])

    dvf = pReg.NiftiImageData3DDisplacement(shift_x, shift_y, shift_z)
    
    return dvf
            

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


def main():

    fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
    input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/pDynamicSimulation/'
    output_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'

    fpath_template_contrast_rawdata = input_fpath_prefix + 'Cube128/CV_nav_cart_128Cube_FLASH_T1.h5'
    fpath_template_acquisition_rawdata = input_fpath_prefix + 'General/meas_MID29_cart_ref_image_FID78804_ismrmrd.h5'

    # configure the simulation
    contrast_rd = pMR.AcquisitionData(fpath_template_contrast_rawdata)
    contrast_rd = pMR.preprocess_acquisition_data(contrast_rd)
    
    acquisition_rd = pMR.AcquisitionData(fpath_template_acquisition_rawdata)
    acquisition_rd = pMR.preprocess_acquisition_data(acquisition_rd)
    
    img_for_contrast = pMR.ImageData()
    img_for_contrast.from_acquisition_data(contrast_rd)

    img_for_acquis = pMR.ImageData()
    img_for_acquis.from_acquisition_data(acquisition_rd)
    img_for_acquis = img_for_acquis.abs()    	

    labels = pReg.NiftiImageData3D( input_fpath_prefix + "Cube128/label_volume.nii"	)
    labels.write(output_fpath_prefix + 'output_example_cartesian_2D_resampling_Labels.nii')

    mvfs = read_motionfields(input_fpath_prefix + 'Cube128/mvf_resp/')
    
    resampler = pReg.NiftyResampler()
    resampler.set_padding_value(0)
    resampler.set_interpolation_type_to_cubic_spline()
    
    resampler.set_floating_image(labels)
    resampler.set_reference_image(labels)

    offset_shift = offset_shift_from_nifti(labels)
    resampler.add_transformation(offset_shift)
    resampler.process()
    
    resampler.clear_transformations()
    resampler.set_floating_image(resampler.get_output())
    resampler.set_reference_image(img_for_acquis)
    resampler.process()

    resampled_img = resampler.get_output()
    resampled_img = resampled_img.abs()
    resampled_img = pReg.NiftiImageData3D(resampled_img)
    
    resampled_img.write(output_fpath_prefix + 'output_example_cartesian_2D_resampling.nii')


    return 1

try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
