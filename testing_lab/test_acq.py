import math
import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'
DATA_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/examples/'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)

from pGadgetron import *

try:
    # acquisitions will be read from this HDF file
    #input_data = MR_Acquisitions(DATA_PATH + 'test_2D_2x.h5')
    input_data = MR_Acquisitions(DATA_PATH + 'testdata.h5')

    processed_data = MR_remove_x_oversampling(input_data)

    na = input_data.number()
    print('%d acquisitions found' % na)

    print('sorting acquisitions...')
    input_data.sort()

    data = abs(input_data.slice_as_array(0))
    pdata = abs(processed_data.slice_as_array(0))
    nx, ny, nc = input_data.slice_dimensions()
    print('acquisition dimensions: %d %d %d' % (nx, ny, nc))
##    shape = data.shape
##    nc = shape[0]
##    ny = shape[1]
##    nx = shape[2]
    for i in range(2):
        pylab.figure(i)
        pylab.imshow(data[i,:,:])
        pylab.figure(i + nc)
        pylab.imshow(pdata[i,:,:])
##        re = abs(re + 1j*im)
##        pylab.imshow(re[i,:,:])
##        minv = numpy.amin(re)
##        maxv = numpy.amax(re)
##        print(minv, maxv)
##        pylab.imshow(re[i,:,:])
##        pylab.figure(i + nc)
##        pylab.imshow(im[i,:,:])
        pylab.show()

    while True:
        s = str(input('enter acquisition number: '))
        if len(s) < 1:
            break
        a = int(s)
        if a < 0 or a >= na:
            break
        acq = input_data.acquisition(a)
        print('flags: %d' % acq.flags())
        print('number of samples: %d' % acq.number_of_samples())
        print('active_channels: %d' % acq.active_channels())
        print('trajectory_dimensions: %d' % acq.trajectory_dimensions())
        print('kspace_encode_step_1: %d' % acq.idx_kspace_encode_step_1())
        print('repetition: %d' % acq.idx_repetition())
        print('slice: %d' % acq.idx_slice())

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
