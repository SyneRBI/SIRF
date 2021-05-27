'''Simplistic PET-MR demo

Usage:
  simplistic_petmr [--help | options]

Options:
  -f <file>, --file=<file>  raw data file
                            [default: simulated_MR_2D_cartesian.h5]
  --mr_path=<path>    path to MR data files, defaults to data/examples/MR
                      subfolder of SIRF root folder
  --mr_engine=<mr>    MR reconstruction engine [default: Gadgetron]
  --pet_engine=<pet>  PET reconstruction engine [default: STIR]
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2015 - 2019 Rutherford Appleton Laboratory STFC
## Copyright 2015 - 2017 University College London.
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

# import common (engine-independent) utilities
import pUtilities as pUtil

# get command-line arguments
__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

# find MR raw data file
data_file = args['--file']
data_path = args['--mr_path']
if data_path is None:
    data_path = pUtil.examples_data_path('MR')
input_file = pUtil.existing_filepath(data_path, data_file)

# import MR and PET engines
exec('import p' + args['--mr_engine' ] + ' as MR' )
exec('import p' + args['--pet_engine'] + ' as PET')

def main():

# MR
    # specify the MR raw data source
    input_data = MR.AcquisitionData(input_file)
    # pre-process acquisitions
    processed_data = MR.preprocess_acquisition_data(input_data)
    # perform reconstruction
    recon = MR.FullySampledReconstructor()
    recon.set_input(processed_data)
    recon.process()
    complex_image = recon.get_output()
    print(complex_image.norm())
    print(complex_image.dot(complex_image))

# PET
    # convert MR image into PET image
    image = PET.ImageData(complex_image)
    print(image.norm())
    print(image.dot(image))
    # apply filter that zeroes the image outside a cylinder of the same
    # diameter as the image xy-section size
    filter = PET.TruncateToCylinderProcessor()
    filter.set_input(image)
    filter.process()
    processed_image = filter.get_output()
    # shortcuts for the above 3 lines
    # image is intact
##    processed_image = filter.process(image)
    # image is modified
##    filter.apply(image)
    # display image
    pUtil.show_3D_array(image.as_array(), \
                        suptitle = 'MR Image', label = 'slice', show = False)
    pUtil.show_3D_array(processed_image.as_array(), \
                        suptitle = 'PET Processed Image', label = 'slice')

try:
    main()
    print('done')
except MR.error as err:
    print('%s' % err.value)
