'''Simplistic PET-MR demo

Usage:
  simplistic_petmr [--help | options]

Options:
  -f <file>, --file=<file>  raw data file
                            [default: simulated_MR_2D_cartesian.h5]
  --mr_path=<path>    path to MR data files, defaults to data/examples/MR
                      subfolder of SIRF root folder
  --mr_engine=<mr>    MR reconstruction engine [default: Gadgetron]
  --pet_engine=<pet>  reconstruction engine [default: Stir]
'''

import pUtil

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)
data_file = args['--file']
data_path = args['--mr_path']
if data_path is None:
    data_path = pUtil.petmr_data_path('mr')
input_file = pUtil.existing_filepath(data_path, data_file)

exec('import p' + args['--mr_engine' ] + ' as MR' )
exec('import p' + args['--pet_engine'] + ' as PET')

def main():

# MR
    # specify the input data file
    input_data = MR.AcquisitionData(input_file)
    # pre-process acquisitions
    processed_data = MR.preprocess_acquisitions(input_data)
    # perform reconstruction
    recon = MR.SimpleReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_image = recon.get_output()

# PET
    # convert MR image into PET image
    image_arr = abs(complex_image.as_array()) # image as Python array
    image = PET.ImageData()                   # empty PET ImageData object
    image.initialise(image_arr.shape[::-1])   # set image shape
    image.fill(image_arr)                     # fill image with values
    # apply cylindric filter
    filter = PET.CylindricFilter()
    filter.apply(image)
    # display image
    image.show()

try:
    main()
except MR.error as err:
    print('exception occured: %s' % err.value)
    
    
