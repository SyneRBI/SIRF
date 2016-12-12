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
('filename', nargs='?', default = 'simulated_MR_2D_cartesian.h5', \
 help = 'raw data file name (default: simulated_MR_2D_cartesian.h5)')
args = parser.parse_args()                                 

def main():

    # acquisitions will be read from an HDF file args.filename
    input_data = AcquisitionData(args.filename)

    na = input_data.number()
    nc, ny, nx = input_data.slice_dimensions()
    print('%d acquisitions found' % na)

    # copy acquisitions into an array
    input_array = input_data.as_array().transpose((1, 0, 2))
    input_shape = input_array.shape
    print('input data dimensions: %dx%dx%d' % input_shape)
    print('input data slice dimensions: %dx%dx%d' % (nc, ny, nx))

    # pre-process acquisition data
    print('processing acquisitions...')
    processed_data = input_data.process(['RemoveROOversamplingGadget'])

    # copy processed acquisitions into an array
    processed_array = processed_data.as_array().transpose((1, 0, 2))
    processed_shape = processed_array.shape
    print('processed data dimensions: %dx%dx%d' % processed_shape)
    print('processed data slice dimensions: %dx%dx%d'\
          % (processed_data.slice_dimensions()))

    print('sorting acquisitions...')
    input_data.sort()

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
        input_slice = abs(input_array[:, (z - 1)*ny : z*ny, :])
        processed_slice = abs(processed_array[:, (z - 1)*ny : z*ny, :])
        print('Enter coil number to view the acquired data for it')
        print('(a value outside the range [1 : %d] will stop this loop)' % nc)
        while True:
            s = str(input('coil: '))
            if len(s) < 1:
                break
            c = int(s)
            if c < 1 or c > nc:
                break
            cp = c + nc
            pylab.figure(c)
            pylab.title('input data')
            pylab.imshow(input_slice[c - 1, :, :])
            pylab.figure(cp)
            pylab.title('processed data')
            pylab.imshow(processed_slice[c - 1, :, :])
            print('Close Figures %d and %d windows to continue...'% (c, cp))
            pylab.show()

##    # perform reconstruction
##    undersampled = input_data.is_undersampled()
##    if not undersampled:
##        print('---\n reconstructing fully sampled data...')
##        recon = SimpleReconstruction()
##    else:
##        print('---\n reconstructing undersampled data using GRAPPA...')
##        recon = GenericCartesianGRAPPAReconstruction()
##    recon.set_input(processed_data)
##    recon.process()
##    if not undersampled:
##        images = recon.get_output()
##    else:
##        images = recon.get_output('images')
##
##    # show obtained images
##    images.show()

try:
    main()
    print('done')

except error as err:
    # display error information
    print('??? %s' % err.value)
