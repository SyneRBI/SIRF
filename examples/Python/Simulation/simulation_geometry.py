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
import nibabel as nib
import time


fpath_testdata_prefix = '/media/sf_CCPPETMR/TestData/'
input_fpath_prefix = fpath_testdata_prefix + 'Input/xDynamicSimulation/pDynamicSimulation/'
output_fpath_prefix = fpath_testdata_prefix + 'Output/xDynamicSimulation/pDynamicSimulation/'

fname_xml = input_fpath_prefix + 'Slab128/XCAT_TissueParameters_XML.xml'
fname_template_contrast = input_fpath_prefix + 'Slab128/CV_nav_cart_128Slab_FLASH_T1.h5'
fname_template_acquisition = input_fpath_prefix + 'General/meas_MID27_CV_11s_TI2153_a6_2x2x8_TR45_FID33312_defaultorientation.h5'
fname_labels = input_fpath_prefix + 'Cube128/label_volume_rai.nii'


reorient_label_volume = False
if reorient_label_volume:

    img = nib.load(fname_labels)
    data = img.get_fdata()

    data = data[:,:,60:70]

    resolution_mm_per_pixel = np.array([2,2,-2,1])
    offset_mm =(-np.array(data.shape)/2 + 0.5) * resolution_mm_per_pixel[0:3]

    affine = np.diag(resolution_mm_per_pixel)
    affine[:3,3] = offset_mm
    
    #
    img = nib.Nifti1Image(data, affine)
    hdr = img.header
    hdr.set_qform(hdr.get_sform())

    fname_out = '/media/sf_CCPPETMR/labels.nii'
    nib.save(img, fname_out)

    sirf_nii = pReg.NiftiImageData(fname_out)
    sirf_nii.print_header()








def resample_to_destination_geometry():

    acquisition_template = pMR.AcquisitionData(fname_template_acquisition)
    labels = pReg.NiftiImageData3D( input_fpath_prefix + "Cube128/label_volume_rai.nii")

    dst_img = pMR.ImageData()
    dst_img.from_acquisition_data(acquisition_template)
    dst_img = pReg.NiftiImageData3D(dst_img)
    
    resampler = pReg.NiftyResampler()
    resampler.set_interpolation_type_to_nearest_neighbour()

    resampler.set_reference_image(dst_img)
    resampler.set_floating_image(labels)

    translation = np.array([0,0,100], dtype=np.float32)
    euler_angles_deg = np.array([0,0,45], dtype=np.float32)

    offset_trafo = pReg.AffineTransformation(translation, euler_angles_deg)
    resampler.add_transformation(offset_trafo)

    resampler.process()
    output = resampler.get_output()
    output.write("/media/sf_CCPPETMR/tmp_neartestneighbor.nii")

def experiments_simulation_geometry():


    contrast_template = pMR.AcquisitionData(fname_template_contrast)
    acquisition_template = pMR.AcquisitionData(fname_template_acquisition)

    labels = pReg.NiftiImageData3D(fname_labels)
    mrsim = pDS.MRDynamicSimulation(labels, fname_xml)

    mrsim.set_contrast_template_data(contrast_template)
    mrsim.set_acquisition_template_data(contrast_template)
    
    # 
    mrsim.save_parametermap_ground_truth(output_fpath_prefix + "simulation_geometry_contrast_parametermap")
    #
    mrsim.set_acquisition_template_data(acquisition_template)
    mrsim.save_parametermap_ground_truth(output_fpath_prefix + "simulation_geometry_acquisition_parametermap")

    #
    offset_x_mm = 0
    offset_y_mm = 0
    offset_z_mm = -10
    
    rotation_angles_deg = [15,15,0]
    translation = np.array([offset_x_mm, offset_y_mm, offset_z_mm])
    euler_angles_deg = np.array(rotation_angles_deg)

    offset_trafo = pReg.AffineTransformation(translation, euler_angles_deg)
    mrsim.set_offset_trafo(offset_trafo)

    mrsim.save_parametermap_ground_truth(output_fpath_prefix + "simulation_geometry_acquisition_offset_parametermap")
    
def main():
    resample_to_destination_geometry()
    # experiments_simulation_geometry()
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
