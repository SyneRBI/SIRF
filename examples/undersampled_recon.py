'''
Upper-level demo, GRAPPA reconstruction of undersampled data.
'''

import argparse
import os
try:
    import pylab
    HAVE_PYLAB = True
except:
    HAVE_PYLAB = False
import sys

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

parser = argparse.ArgumentParser(description = \
'''
Upper-level demo, GRAPPA reconstruction of undersampled data.
''')
parser.add_argument\
('filename', nargs='?', default = 'simulated_MR_2D_cartesian_Grappa2.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian_Grappa2.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)
    if not input_data.is_undersampled():
        print('this demo needs undersampled raw data')
        return

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoAdjustROGadget', \
         'RemoveROOversamplingGadget']
    print('---\n pre-processing acquisitions...')
    preprocessed_data = input_data.process(prep_gadgets)

    # perform reconstruction
    recon = GenericCartesianGRAPPAReconstruction()
    # for undersampled acquisition data GRAPPA will compute Gfactor images
    # in addition to reconstructed ones
    recon.compute_gfactors(True)
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    image = recon.get_output('image')
    gfactor = recon.get_output('gfactor')
    data = abs(image.as_array())
    gdata = abs(gfactor.as_array())

    nz = data.shape[0]
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
        pylab.imshow(data[z - 1,:,:])
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
