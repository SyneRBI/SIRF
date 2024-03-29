'''FBP2D reconstruction demo

Usage:
  fbp2d_reconstruction [--help | options]

Options:
  -d <file>, --file=<file>    raw data file
                              [default: simulated_data.hs]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/PET
                              subfolder of SIRF root folder
  -e <engn>, --engine=<engn>  reconstruction engine [default: STIR]
  --non-interactive           do not show plots
'''

## SyneRBI Synergistic Image Reconstruction Framework (SIRF)
## Copyright 2018 - 2020 Rutherford Appleton Laboratory STFC
## Copyright 2018 University College London.
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

from sirf.Utilities import error, examples_data_path, existing_filepath

# import engine module
import importlib
engine = args['--engine']
pet = importlib.import_module('sirf.' + engine)


# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('PET')
raw_data_file = existing_filepath(data_path, data_file)
show_plot = not args['--non-interactive']


def main():

    # no info printing from the engine, warnings and errors sent to stdout
    _ = pet.MessageRedirector()

    # PET acquisition data to be read from the file specified by --file option
    print('raw data: %s' % raw_data_file)
    acq_data = pet.AcquisitionData(raw_data_file)

    # create reconstructor object
    recon = pet.FBP2DReconstructor()
    # specify the acquisition data
    recon.set_input(acq_data)

    # reconstruct with default settings
    recon.process()
    image = recon.get_output()
    nz, ny, nx = image.dimensions()
    print('--------\n dimensions: nx = %d, ny = %d, nz = %d' % (nx, ny, nz))
    z = int(nz*2/3)
    if show_plot:
        image.show(z)

    # change image size
    recon.set_output_image_size_xy(nx*2)
    recon.process()
    image = recon.get_output()
    nz, ny, nx = image.dimensions()
    print('--------\n dimensions: nx = %d, ny = %d, nz = %d' % (nx, ny, nz))
    if show_plot:
        image.show(z)

    # zoom in
    zoom = 2.5
    recon.set_zoom(zoom)
    recon.process()
    image = recon.get_output()
    print('--------\n zoom %f' % zoom)
    if show_plot:
        image.show(z)

    # use a Hann filter
    alpha = 0.5
    recon.set_alpha_cosine_window(alpha)
    recon.process()
    image = recon.get_output()
    print('--------\n alpha %f' % alpha)
    if show_plot:
        image.show(z)

    # a Hann filter with lower cut-off (0.5 is no cut-off)
    fc = 0.2
    recon.set_frequency_cut_off(fc)
    recon.process()
    image = recon.get_output()
    print('--------\n frequency cut-off %f' % fc)
    if show_plot:
        image.show(z)

    # alternative way to set the output image parameters (via image template)
    image_tmpl = acq_data.create_uniform_image() # image template
    recon.set_up(image_tmpl) # use image template to create the output image
    recon.process()
    image = recon.get_output()
    print('--------\n alternative setup')
    if show_plot:
        image.show(z)


# if anything goes wrong, an exception will be thrown 
# (cf. Error Handling section in the spec)
try:
    main()
    print('\n=== done with %s' % __file__)

except error as err:
    # display error information
    print('%s' % err.value)
