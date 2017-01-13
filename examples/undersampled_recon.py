'''
Upper-level demo, GRAPPA reconstruction of undersampled data.

Usage:
  undersampled_recon.py [--help | options]

Options:
  -f <file>, --file=<file>    raw data file
                              [default: simulated_MR_2D_cartesian_Grappa2.h5]
  -e <engn>, --engine=<engn>  reconstruction engine [default: Gadgetron]
  -p <path>, --path=<path>    sub-path to engine module
                              [default: /xGadgetron/pGadgetron]
'''

__version__ = '0.1.0'
from docopt import docopt
args = docopt(__doc__, version=__version__)

import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys

sys.path.append(os.environ.get('SRC_PATH') + args['--path'])

exec('from p' + args['--engine'] + ' import *')

def main():

    # acquisitions will be read from an HDF file
    input_data = AcquisitionData(args['--file'])
    if not input_data.is_undersampled():
        print('this demo needs undersampled raw data')
        return

    # pre-process acquisitions
    print('---\n pre-processing acquisitions...')
    preprocessed_data = preprocess_acquisitions(input_data)

    # perform reconstruction
    recon = GenericCartesianGRAPPAReconstruction()
    # for undersampled acquisition data GRAPPA will compute Gfactor images
    # in addition to reconstructed ones
    recon.compute_gfactors(True)
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    image = recon.get_output('image')
    gfact = recon.get_output('gfactor')
    idata = abs(image.as_array())
    gdata = abs(gfact.as_array())

    nz = idata.shape[0]
    # plot image and gfactor slices
    while HAVE_PYLAB:
        print('---\n Enter the slice number to view it.')
        print(' A value outside the range [1 : %d] will stop this loop.'% nz)
        s = str(input('slice: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 1 or z > nz:
            break
        pylab.figure(z)
        pylab.title('image')
        pylab.imshow(idata[z - 1,:,:])
        print('Close Figure %d window to continue...' % z)
        pylab.figure(z + nz)
        pylab.title('G factor')
        pylab.imshow(gdata[z - 1,:,:])
        print('Close Figure %d window to continue...' % (z + nz))
        pylab.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
