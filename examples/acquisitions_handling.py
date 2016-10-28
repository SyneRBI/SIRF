'''
Upper-level interface demo that illustrates acquisitions
pre-processing, sorting and plotting.
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
Upper-level interface demo that illustrates acquisitions
pre-processing, sorting and plotting.
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = MR_Acquisitions(args.filename)

    # pre-process acquisition data
    print('processing acquisitions...')
    prep_gadgets = ['NoiseAdjustGadget', 'AsymmetricEchoGadget', \
         'RemoveROOversamplingGadget']
    processed_data = input_data.process(prep_gadgets)

    na = input_data.number()
    print('%d acquisitions found' % na)

    print('sorting acquisitions...')
    input_data.sort()

    nx, ny, nc = input_data.slice_dimensions()
    nz = na//ny

    while HAVE_PYLAB:
        print('---\n Enter the slice number to view it.')
        print(' A value outside the range [1 : %d] will stop this loop.'% nz)
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 1 or z > nz:
            break
        data = abs(input_data.slice_as_array(z - 1))
        pdata = abs(processed_data.slice_as_array(z - 1))
        print('Enter coil number to view the acquired data for it')
        print('(a value outside the range [1 : %d] will stop this loop)' % nc)
        while True:
            s = str(input('coil: '))
            if len(s) < 1:
                break
            c = int(s)
            if c < 1 or c > nc:
                break
            pylab.figure(c)
            pylab.title('input data')
            pylab.imshow(data[c - 1, :, :])
            pylab.figure(c + nc)
            pylab.title('processed data')
            pylab.imshow(pdata[c - 1, :, :])
            print('Close Figures %d and %d windows to continue...'% (c, c + nc))
            pylab.show()

    # perform reconstruction
    recon = MR_BasicReconstruction()
    recon.set_input(processed_data)
    recon.process()
    complex_images = recon.get_output()

    # extract real images from complex
    images = complex_images.real()
    # show obtained images
    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
