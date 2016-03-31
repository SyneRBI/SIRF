import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *
from pGadgets import *

try:
    # acquisitions will be read from this HDF file
    input_data = ISMRMRDAcquisitions('testdata.h5')
    #input_data = ISMRMRDAcquisitions('opismrmrd.h5')

    na = input_data.number()
    print('%d acquisitions found' % na)
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
