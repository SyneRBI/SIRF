'''Simplistic PET-MR demo

Usage:
  using_acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian.h5]
  -p <path>, --path=<path>    path to data files, defaults to data/examples/MR
                              subfolder of SIRF root folder
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import pGadgetron as MR
import pStir as PET

def main():

    # locate the input data file
    data_path = args['--path']
    if data_path is None:
        data_path = MR.mr_data_path()
    input_file = MR.existing_filepath(data_path, args['--file'])

    # acquisitions will be read from an HDF file input_file
    input_data = MR.AcquisitionData(input_file)

    print('---\n acquisition data norm: %e' % input_data.norm())

    # pre-process acquisitions
    print('---\n processing acquisitions...')
    processed_data = MR.preprocess_acquisitions(input_data)

    print('---\n processed acquisition data norm: %e' % processed_data.norm())

    # perform reconstruction
    recon = MR.SimpleReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    image_arr = abs(complex_images.as_array())
    nz, ny, nx = image_arr.shape
    image = PET.ImageData()
    image_size = (nx, ny, nz)
    voxel_size = (3, 3, 3.375)
    image.initialise(image_size, voxel_size)
    image.fill(image_arr)
    image.show()

try:
    main()
except MR.error as err:
    print('exception occured: %s' % err.value)
    
    
