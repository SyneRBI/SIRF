'''
Upper-level demo that illustrates the computation of coil sensitivity maps
and applying projection from the image space into acquisition space and back
defined by the aquisition model.

Usage:
  acquisition_model.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
  -o <file>, --output=<file>  output file for simulated data
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2017 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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

# import engine module
exec('from p' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
output_file = args['--output']

def main():

    # locate the k-space raw data file
    input_file = existing_filepath(data_path, data_file)

    # acquisition data will be read from an HDF file input_file
    acq_data = AcquisitionData(input_file)
    print('---\n acquisition data norm: %e' % acq_data.norm())

    # pre-process acquisition data
    print('---\n pre-processing acquisition data...')
    processed_data = preprocess_acquisition_data(acq_data)
    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction to obtain a meaningful ImageData object
    # (cannot be obtained in any other way at present)
#    if processed_data.is_undersampled():
    if acq_data.is_undersampled():
        recon = CartesianGRAPPAReconstructor();
        recon.compute_gfactors(False)
    else:
        recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    print('---\n reconstructing...')
    recon.process()
    reconstructed_images = recon.get_output()
    r_norm = reconstructed_images.norm()
    print('---\n reconstructed images norm: %e' % r_norm)

    for i in range(min(8, reconstructed_images.number())):
        reconstructed_image = reconstructed_images.image(i)
        print('--- image %d' % i)
        for p in [ \
            'version', 'flags', 'data_type', 'channels', \
            'slice', 'repetition', \
            'image_type', 'image_index', 'image_series_index' \
            ]:
            form = p + ' %d'
            print(form % reconstructed_image.info(p))
        print('matrix size:'),
        print(reconstructed_image.matrix_size())
        print('patient_table_position:'),
        print(reconstructed_image.patient_table_position())

    ind = reconstructed_images.get_info('image_index')
    print('image indices:')
    print(ind)
    ind = reconstructed_images.get_info('slice')
    print('image slices:')
    print(ind)
    ptp = reconstructed_images.get_info('patient_table_position')
    print('patient table positions:')
    print(ptp)
    
    # sort processed acquisition data;
    # sorting currently performed with respect to (in this order):
    #    - repetition
    #    - slice
    #    - kspace encode step 1
    print('---\n sorting acquisition data...')
    processed_data.sort()
    if output_file is not None:
        processed_data.write(output_file)

    # compute coil sensitivity maps
    print('---\n computing coil sensitivity maps...')
    csms = CoilSensitivityData()
    csms.calculate(processed_data)
    # alternatively, coil sensitivity maps can be computed from
    # CoilImageData - see coil_sensitivity_maps.py

    # create acquisition model based on the acquisition parameters
    # stored in processed_data and image parameters stored in reconstructed_images
##    acq_model = AcquisitionModel(processed_data, reconstructed_images)
    acq_model = AcquisitionModel()
    acq_model.set_up(processed_data, reconstructed_images)
    acq_model.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce simulated
    # acquisition data
    print('---\n forward-projecting...')
    simulated_acq_data = acq_model.forward(reconstructed_images)
    print('---\n reconstructed images forward projection norm %e'\
          % simulated_acq_data.norm())
    if output_file is not None:
        simulated_acq_data.write(output_file)

    # display simulated acquisition data
    #simulated_acq_data.show(title = 'Simulated acquisition data (magnitude)')

    # backproject simulated acquisition data
    print('---\n backprojecting...')
    backprojected_data = acq_model.backward(simulated_acq_data)
    b_norm = backprojected_data.norm()
    print('norm of backprojected images: %f' % b_norm)
    # display backprojected data
    backprojected_data.show(title = 'Backprojected data (magnitude)')
    # display reconstructed images
    reconstructed_images.show(title = 'Reconstructed images (magnitude)')

    diff = backprojected_data/b_norm - reconstructed_images/r_norm
    print('norm of backprojected - reconstructed images: %f' % diff.norm())
    # testing fill
    reconstructed_images.fill(backprojected_data)
    diff = (backprojected_data - reconstructed_images).norm()
    if diff > 0:
        print('fill error: %f' % diff)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
