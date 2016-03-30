import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *

try:
    csms = MRCoilSensitivityMaps('csm_opismrmrd.h5')
    nz = csms.number()
    print('%d slices' % nz)

    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = csms.csm_as_array(z)
        shape = data.shape
        nc = shape[0]
        for i in range(nc):
            pylab.figure(z*nc + i + 1)
            pylab.imshow(data[i,0,:,:])
            pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
