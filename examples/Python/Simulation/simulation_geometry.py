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
fname_template_acquisition = input_fpath_prefix + 'General/meas_MID33_rad_2d_gc_FID78808_ismrmrd_defaultorient.h5'
fname_labels = input_fpath_prefix + 'Slab128/label_volume_rai.nii'

def reorient_label_volume():

    img = nib.load(fname_labels)
    data = img.get_fdata()

    data = data[:,:,60:70]

    affine = np.diag([2,2,-2,1])
    affine[:,3] = [-127, -127, 9 , 1]

    img = nib.Nifti1Image(data, affine)
    fname_out = input_fpath_prefix + 'Slab128/label_volume_rai.nii'

    nib.save(img, fname_out)

# reorient_label_volume()


acquisition_template = pMR.AcquisitionData(fname_template_acquisition)
contrast_template = pMR.AcquisitionData(fname_template_contrast)

def experiments_simulation_geometry():

    labels = pReg.NiftiImageData3D(fname_labels)
    mrsim = pDS.MRDynamicSimulation(labels, fname_xml)

    mrsim.set_contrast_template_data(contrast_template)
    mrsim.set_acquisition_template_data(contrast_template)
    
    # 
    mrsim.save_parametermap_ground_truth(output_fpath_prefix + "simulation_geometry_contrast_parametermap_")
    #
    mrsim.set_acquisition_template_data(acquisition_template)
    mrsim.save_parametermap_ground_truth(output_fpath_prefix + "simulation_geometry_acquisition_parametermap_")

    #
    offset_x_mm = 0
    offset_y_mm = 0
    offset_z_mm = -10
    
    rotation_angles_deg = [15,15,0]
    translation = np.array([offset_x_mm, offset_y_mm, offset_z_mm])
    euler_angles_deg = np.array(rotation_angles_deg)

    offset_trafo = pReg.AffineTransformation(translation, euler_angles_deg)
    mrsim.set_offset_trafo(offset_trafo)

    mrsim.save_parametermap_ground_truth(output_fpath_prefix + "simulation_geometry_acquisition_offset_parametermap_")
    
def main():
    # print_header_infos()
    experiments_simulation_geometry()
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('??? %s' % err.value)
    exit(1)
