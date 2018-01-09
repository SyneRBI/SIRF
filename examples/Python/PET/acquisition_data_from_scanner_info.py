'''Acquisition data from scanner info demo.

Usage:
  acquisition_data_from_scanner_info [--help | options]

Options:
  -e <engn>, --engine=<engn>   reconstruction engine [default: STIR]
  -s <stsc>, --storage=<stsc>  acquisition data storage scheme [default: file]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

from pUtilities import show_2D_array

# import engine module
exec('from p' + args['--engine'] + ' import *')

storage = args['--storage']

def main():

    AcquisitionData.set_storage_scheme(storage)

    acq_data = AcquisitionData('Siemens_mMR')

    # copy the acquisition data into a Python array
    acq_array = acq_data.as_array()
    acq_dim = acq_array.shape
    print(acq_dim)
    z = acq_dim[0]//2

    show_2D_array('Acquisition data', acq_array[z,:,:])

try:
    main()
    print('done')
except error as err:
    print('%s' % err.value)
