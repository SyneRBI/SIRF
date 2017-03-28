'''
Upper-level demo that illustrates the computation of coil sensitivity maps
and applying projection from the image space into acquisition space and back
defined by the aquisition model.

Usage:
  simple_simulation.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
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

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = mr_data_path()
    input_file = existing_filepath(data_path, args['--file'])

    # acquisitions will be read from an HDF file input_file
    input_data = AcquisitionData(input_file)

    print('---\n acquisition data norm: %e' % input_data.norm())

    # pre-process acquisition data
    print('---\n pre-processing acquisition data...')
    processed_data = preprocess_acquisition_data(input_data)

    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction to obtain a meaningful ImageData object
    # (cannot be obtained in any other way at present)
    recon = FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    print('---\n reconstructed images norm: %e' % complex_images.norm())

    csms = CoilSensitivityData()

    print('---\n sorting acquisitions...')
    processed_data.sort()
    print('---\n computing coil sensitivity maps...')
    csms.calculate(processed_data)
    # alternatively, coil sensitivity maps can be computed from
    # CoilImageData - see coil_sensitivity_maps.py

    # create acquisition model based on the acquisition parameters
    # stored in input_data and image parameters stored in complex_images
    acq_model = AcquisitionModel(processed_data, complex_images)

    acq_model.set_coil_sensitivity_maps(csms)

    # use the acquisition model (forward projection) to produce simulated
    # acquisition data
    simulated_acq_data = acq_model.forward(complex_images)

    print('---\n reconstructed images forward projection norm %e'\
          % simulated_acq_data.norm())

    # get data as a Python ndarray
    simulated_acq_array = simulated_acq_data.as_array();

    # TODO display a slice etc

    # backproject simulated acquisition data
    backprojected_data = acq_model.backward(simulated_acq_data)

    # show backprojected images
    backprojected_data.show()


try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
