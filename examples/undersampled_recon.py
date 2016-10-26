'''
Upper-level demo, GRAPPA reconstruction of undersampled data.
'''

import argparse
import os
import pylab
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
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)
    if not input_data.is_undersampled():
        print('this demo needs undersampled raw data')
        return

    # pre-process acquisitions
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoGadget', \
         'RemoveROOversamplingGadget']
    print('---\n pre-processing acquisitions...')
    preprocessed_data = input_data.process(prep_gadgets)

    # perform reconstruction
    recon = MR_BasicGRAPPAReconstruction()
    # for undersampled acquisition data GRAPPA will compute Gfactor images
    # in addition to reconstructed ones
    recon.compute_gfactors(True)
    recon.set_input(preprocessed_data)
    print('---\n reconstructing...')
    recon.process()
    complex_images, complex_gfactors = recon.get_output()

    # get real-valued reconstructed images and gfactors
    print('---\n processing images...')
    images = complex_images.real()
    gfactors = complex_gfactors.real()

    nz = images.number()
    print('%d images reconstructed.' % nz)

    # plot images and gfactors
    print('Enter the number of the slice to view it')
    print('(a value outside the range [1 : %d] will stop this loop)'% nz)
    while True:
        s = str(input('slice: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 1 or z > nz:
            break
        data = images.image_as_array(z - 1)
        gdata = gfactors.image_as_array(z - 1)
        pylab.figure(z)
        pylab.title('image')
        pylab.imshow(data[0,0,:,:])
        print('Close Figure %d window to continue...' % z)
        pylab.figure(z + nz)
        pylab.title('G factor')
        pylab.imshow(gdata[0,0,:,:])
        print('Close Figure %d window to continue...' % (z + nz))
        pylab.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
