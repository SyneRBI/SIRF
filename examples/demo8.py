'''
Upper level interface demo that illustrates acquisitions
pre-processing, sorting and plotting.
'''

import math
import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    input_data = MR_Acquisitions('testdata.h5')

    # pre-process acquisition data
    print('processing acquisitions...')
    processed_data = MR_remove_x_oversampling(input_data)

    na = input_data.number()
    print('%d acquisitions found' % na)

    print('sorting acquisitions...')
    input_data.sort()

    nx, ny, nc = input_data.dimensions()
    nz = na//ny

    print('Enter z-coordinate of the slice to view the acquired data for it')
    print('(a value outside the range [0 : %d) will stop this loop)' % nz)
    while True:
        s = str(input('z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = abs(input_data.slice_as_array(z))
        pdata = abs(processed_data.slice_as_array(z))
        print('Enter coil number to view the acquired data for it')
        print('(a value outside the range [0 : %d) will stop this loop)' % nc)
        while True:
            s = str(input('coil: '))
            if len(s) < 1:
                break
            c = int(s)
            if c < 0 or c >= nc:
                break
            pylab.figure(c)
            pylab.title('oversampled data')
            pylab.imshow(data[c,:,:])
            pylab.figure(c + nc)
            pylab.title('de-oversampled data')
            pylab.imshow(pdata[c,:,:])
            print('Close Figures %d and %d windows to continue...'% (c, c + nc))
            pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
