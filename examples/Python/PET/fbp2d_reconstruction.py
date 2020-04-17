'''FBP2D reconstruction demo

Usage:
  fbp2d_reconstruction [--help | options]

Options:
  -d <file>, --file=<file>    raw data file
                              [default: my_forward_projection.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
'''

## CCP PETMR Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 Rutherford Appleton Laboratory STFC
## Copyright 2018 University College London.
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
exec('from sirf.' + args['--engine'] + ' import *')

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)

def main():

    # no info printing from the engine, warnings and errors sent to stdout
    msg_red = MessageRedirector()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = AcquisitionData(raw_data_file)

    # create reconstructor object
    recon = FBP2DReconstructor()
    # specify the acquisition data
    recon.set_input(acq_data)

    # reconstruct with default settings
    recon.process()
    image = recon.get_output()
    image_array = image.as_array()
    z = int(image_array.shape[0]*2/3)
    print('--------\n xy-size %d' % image_array.shape[1])
    image.show(z)

    # change image size
    recon.set_output_image_size_xy(image_array.shape[1]*2)
    recon.process()
    image = recon.get_output()
    image_array = image.as_array()
    print('--------\n xy-size %d' % image_array.shape[1])
    image.show(z)

    # zoom in
    zoom = 2.5
    recon.set_zoom(zoom)
    recon.process()
    image = recon.get_output()
    print('--------\n zoom %f' % zoom)
    image.show(z)

    # use a Hann filter
    alpha = 0.5
    recon.set_alpha_cosine_window(alpha)
    recon.process()
    image = recon.get_output()
    print('--------\n alpha %f' % alpha)
    image.show(z)

    # a Hann filter with lower cut-off (0.5 is no cut-off)
    fc = 0.2
    recon.set_frequency_cut_off(fc)
    recon.process()
    image = recon.get_output()
    print('--------\n frequency cut-off %f' % fc)
    image.show(z)

    # alternative way to set the output image parameters (via image template)
    image1 = acq_data.create_uniform_image() # image template
    recon.set_up(image1) # use image template to create the output image
    recon.process()
    image = recon.get_output()
    print('--------\n alternative setup')
    image.show(z)

# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('--------\n done')
except error as err:
    # display error information
    print('%s' % err.value)
