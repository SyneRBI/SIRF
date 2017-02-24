'''Simplistic PET-MR demo

Usage:
  using_acquisition_model [--help | options]

Options:
  -f <file>, --file=<file>   raw data file
                                  [default: simulated_MR_2D_cartesian.h5]
  --mr_path=<path>      path to MR data files, defaults to data/examples/MR
                                  subfolder of SIRF root folder
  --mr_engine=<mr>     MR reconstruction engine [default: Gadgetron]
  --pet_engine=<pet>    reconstruction engine [default: Stir]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

exec('import p' + args['--mr_engine'] + ' as MR')
exec('import p' + args['--pet_engine'] + ' as PET')

def main():

    # locate the input data file
    data_path = args['--mr_path']
    if data_path is None:
        data_path = MR.mr_data_path()
    input_file = MR.existing_filepath(data_path, args['--file'])

    # acquisitions will be read from an HDF file input_file
    input_data = MR.AcquisitionData(input_file)

    # pre-process acquisitions
    processed_data = MR.preprocess_acquisitions(input_data)

    # perform reconstruction
    recon = MR.SimpleReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_image = recon.get_output()

    # convert MR image into PET image
    image = PET.ImageData()
    image_arr = abs(complex_image.as_array())
    nz, ny, nx = image_arr.shape
    image_size = (nx, ny, nz)
    voxel_size = (3, 3, 3.375)
    image.initialise(image_size, voxel_size)
    image.fill(image_arr)

    # apply cylindric filter
    filter = PET.CylindricFilter()
    filter.apply(image)

    image.show()

try:
    main()
except MR.error as err:
    print('exception occured: %s' % err.value)
    
    
