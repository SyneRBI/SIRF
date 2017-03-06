'''Simplistic PET-MR demo

Usage:
  using_acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
'''

import pUtil
import pGadgetron as MR # define MR engine
import pStir as PET     # define PET engine

# process command-line options
__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)
data_file = args['--file']
data_path = args['--path']
if data_path is None:
    data_path = pUtil.petmr_data_path('mr')
input_file = pUtil.existing_filepath(data_path, data_file)


def main():

# MR
    # specify the source of raw MR data
    input_data = MR.AcquisitionData(input_file)
    # pre-process MR data
    processed_data = MR.preprocess_acquisitions(input_data)
    # perform MR reconstruction
    recon = MR.SimpleReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_image = recon.get_output()

# PET
    # convert MR image into PET image
    image_arr = abs(complex_image.as_array()) # image as Python array
    image = PET.ImageData() # empty PET ImageData object
    image.initialise(image_arr.shape[::-1]) # set image shape
    image.fill(image_arr) # fill image with values from image_array
    # apply cylindric filter
    filter = PET.CylindricFilter()
    filter.apply(image)
    # display image
    image.show()

try:
    main()
except MR.error as err:
    print('exception occured: %s' % err.value)
    
    
