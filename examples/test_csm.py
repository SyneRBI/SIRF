import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *

try:
    #input_data = ISMRMRDAcquisitions('testdata.h5')
    #input_data = ISMRMRDAcquisitions('opismrmrd.h5')
    input_data = ISMRMRDAcquisitions('nn_no.h5')

    print('ordering acquisitions...')
    input_data.order()

    print('computing sensitivity maps...')
    csms = MRCoilSensitivityMaps()
    csms.compute(input_data)

    nz = csms.number()
    print('%d slices' % nz)

    maxv = 0
    for z in range(nz):
        data = csms.csm_as_array(z)
        minvz = numpy.amin(data)
        maxvz = numpy.amax(data)
        if z == 0:
            minv = minvz
        else:
            minv = min(minvz, minv)
        maxv = max(maxvz, maxv)
    print(minv, maxv)

    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        data = csms.csm_as_array(z)/maxv
        shape = data.shape
        nc = shape[0]
        for i in range(nc):
            pylab.figure(z*nc + i + 1)
            pylab.imshow(data[i,0,:,:], vmin = 0, vmax = 1)
            #pylab.imshow(data[i,0,:,:])
            pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
