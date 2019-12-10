'''
Medium-level interface demo that illustrates 2D Cartesian MR image 
reconstruction using Gadgetron by directly creating and running a chain of 
gadgets.

Usage:
  fully_sampled_recon_single_chain.py [--help | options]

Options:
  -f <file>, --file=<file>           raw data file
                                     [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>           path to data files, defaults to data/examples/MR
                                     subfolder of SIRF root folder
  -o <file>, --output=<file>         images output file
  -a <string>, --algorithm=<string>  algorithm to use ('SimpleReconGadget', 'GenericReconCartesianFFTGadget') [default: SimpleReconGadget]
  --type_to_save=<string>            type to save ('mag', 'imag', 'all') [default: all]
  --show                             show plots
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

import time
import sys

# import SIRF utilities
from sirf.Utilities import examples_data_path, existing_filepath, error
# import MR engine types
from sirf.Gadgetron import AcquisitionData, Reconstructor

# process command-line options
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = examples_data_path('MR')
output_file = args['--output']

type_to_save = args['--type_to_save']
show_plot = False
if args['--show']:
    show_plot = True

algorithm = args['--algorithm']

def main():

    # locate the input data
    input_file = existing_filepath(data_path, data_file)
    acq_data = AcquisitionData(input_file)
    
    if algorithm == 'SimpleReconGadget':
        extra_gadgets = [algorithm]
    else:
        extra_gadgets = [algorithm, 'GenericReconFieldOfViewAdjustmentGadget']
    
    # create reconstruction object
    # Rather than using a predefined image reconstruction object, here a new 
    # image reconstruction object is created by concatinating multiple gadgets 
    # (for more information on Gadgetron and its gadgets please see: 
    # https://github.com/gadgetron/.).
    # Parameters for individual gadgets can be defined either during the 
    # creation of the reconstruction object:
    #   e.g. AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)
    # or by giving a gadget a label (cf. label ex: for the last gadget)
    # and using set_gadget_property(label, propery, value).
    # The gadgets will be concatenated and will be executed as soon as 
    # process() is called.
    recon_gadgets = ['NoiseAdjustGadget',
        'AsymmetricEchoAdjustROGadget',
        'RemoveROOversamplingGadget',
        'AcquisitionAccumulateTriggerGadget(trigger_dimension=repetition)',
        'BucketToBufferGadget(split_slices=true, verbose=false)'] \
        + extra_gadgets + \
        ['ImageArraySplitGadget', 
        'ex:ExtractGadget'
        ]

    recon = Reconstructor(recon_gadgets)

    # ExtractGadget defines which type of image should be returned:
    # none      0
    # magnitude 1
    # real      2
    # imag      4
    # phase     8
    # in this example '5' returns both magnitude and imaginary part
##    recon.set_gadget_property('ex', 'extract_mask', 5)
    # === THE ABOVE IS OBSOLETE, NOW SHOULD USE ===>
    if type_to_save=='mag' or type_to_save=='all':
        recon.set_gadget_property('ex', 'extract_magnitude', True)
    if type_to_save=='imag' or type_to_save=='all':
        recon.set_gadget_property('ex', 'extract_imag', True)
    
    # provide raw k-space data as input
    recon.set_input(acq_data)

    # optionally set Gadgetron server host and port
    recon.set_host('localhost')
    # On VM you can try a port other than the default 9002, e.g. 9003, by taking
    # the following steps:
    # 1) in ~/devel/install/share/gadgetron/config/gadgetron.xml replace
    #    <port>9002</port> with <port>9003</port>
    # 2) go to Settings->Network->Advanced->Port Forwarding and add new rule
    #    (click on green + in the upper right corner) with Host and Guest ports
    #    set to 9003
    # 3) uncomment the next line
    #recon.set_port('9003')
    # Note: each gadget chain can run on a different VM - to try, start two VMs
    # and do the above steps 1 and 2 on one of them, then add
    # recon.set_port('9003') before recon.process in grappa_detail.py
    # (where preprocessing will still run on default port 9002). 

    # perform reconstruction
    recon.process()
    
    # retrieve reconstructed image data
    image_data = recon.get_output()

    # show reconstructed image data
    if show_plot:
        for im in range(image_data.number()):
            image = image_data.image(im)
            # image types   series
            # magnitude 1       0
            # phase     2    3000
            # real      3    1000
            # imag      4    2000
            im_type = image.image_type()
            im_series = image.image_series_index()
            print('image: %d, type: %d, series: %d' % (im, im_type, im_series))
        image_data.show(title = 'Images magnitude and imaginary part')

    if output_file is not None:
        # write images to a new group in args.output
        # named after the current date and time
        time_str = time.asctime()
        print('writing to %s' % output_file)
        image_data.write(output_file) #, time_str)

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
    sys.exit(1)
