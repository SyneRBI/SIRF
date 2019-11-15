'''
Medium-level demo demonstrating how 2D coil sensitivity maps can be obtained 
from a multi-coil 2D Cartesian MR acquisition

Usage:
  coil_sensitivity_maps.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -i <iter>, --iter=<iter>    number of smoothing iterations [default: 10]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF).
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC.
## Copyright 2015 - 2017 University College London.
## Copyright 2015 - 2017 Physikalisch-Technische Bundesanstalt.
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

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from pUtilities import *

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
nit = int(args['--iter'])

def main():

    # 1. Prepare data:

    # locate the input data file
    input_file = existing_filepath(data_path, data_file)
    #
    # acquisition data will be read from an HDF file input_file
    acq_data = AcquisitionData(input_file)
    #
    # pre-process acquisition data
    processed_data = preprocess_acquisition_data(acq_data)
    #
    # sort k-space data into a 2D Cartesian matrix for each coil
    processed_data.sort()
    
    # 2. Calculate coil sensitivity maps directly from the raw k-space data:

    # create coil sensitivity object
    CSMs = CoilSensitivityData()
    #
    # set number of smoothing iterations to suppress noise
    CSMs.smoothness = nit
    #
    # calculate coil sensitivity maps directly from the raw k-space data by the
    # Square-Root-of-the-Sum-of-Squares over all coils (SRSS)
    CSMs.calculate(processed_data)
    #
    # display coil sensitivity maps
    csms_array = CSMs.as_array()
    nz = csms_array.shape[1]
    title = 'SRSS from raw data (magnitude)'
    show_3D_array(abs(csms_array[:, nz//2, :, :]), suptitle=title, \
                  xlabel='samples', ylabel='readouts', label='coil', show=False)

    # 3. Now compute coil sensitivity maps from coil images in order to compare
    # SSRS and Inati methods:

    # create object containing images for each coil
    CIs = CoilImageData()
    #
    # calculate coil images from raw data
    CIs.calculate(processed_data)
    #
    # create coil sensitivity object
    CSMs = CoilSensitivityData()
    #
    # calculate coil sensitivity maps by dividing each coil image data by the
    # Square-Root-of-the-Sum-of-Squares over all coils (SRSS);
    # (niter = nit) sets the number of smoothing iterations applied
    # to the image data prior to the calculation of the coil sensitivity maps
    CSMs.calculate(CIs, method='SRSS(niter = %d)' % nit)
    #
    # display coil sensitivity maps (must be identical to previously computed)
    csms_array = CSMs.as_array()
    nz = csms_array.shape[1]
    title = 'SRSS from coil images (magnitude)'
    show_3D_array(abs(csms_array[:, nz//2, :, :]), suptitle=title, \
                  xlabel='samples', ylabel='readouts', label='coil', \
                  show=False)

    try:
        from ismrmrdtools import coils
        # calculate coil sensitivity maps using an approach suggested by 
        #   Inati SJ, Hansen MS, Kellman P.
        #   A solution to the phase problem in adaptive coil combination.
        #   In: ISMRM proceeding; April; Salt Lake City, Utah, USA; 2013. 2672.  
        # for more details please see 
        # gadgetron/toolboxes/mri_core/mri_core_coil_map_estimation.h  
        CSMs = CoilSensitivityData()
        CSMs.calculate(CIs, method='Inati()')
        csms_array = CSMs.as_array()
        #
        # display coil sensitivity maps
        title = 'Inati (magnitude)'
        show_3D_array(abs(csms_array[:, nz//2, :, :]), suptitle=title, \
                      xlabel='samples', ylabel='readouts', label='coil')
    except:
        print('ismrmrd-python-tools not found, skipping Inati method')

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
