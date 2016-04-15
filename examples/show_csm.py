import math
import pylab
import sys
import time

sys.path.append('../../build/xGadgetron')
sys.path.append('../pGadgetron')

from pGadgetron import *

try:
    csms = MR_CoilSensitivityMaps()
    csm_file = str(input('csm file: '))
    print('reading sensitivity maps...')
    csms.read(csm_file)

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
        #shape = data.shape
        re, im = csms.csm_as_arrays(z)/maxv        
        shape = re.shape
        nc = shape[0]
        ny = shape[1]
        nx = shape[2]
        for i in range(nc):
            for iy in range(ny):
                for ix in range(nx):
                    im[i,0,iy,ix] = math.atan2(im[i,0,iy,ix], re[i,0,iy,ix])
            pylab.figure(z*nc + i + 1)
            pylab.imshow(data[i,0,:,:], vmin = 0, vmax = 1)
##            pylab.figure((z + 1)*nc + i + 1)
##            pylab.imshow(re[i,0,:,:], vmin = -1, vmax = 1)
##            pylab.colorbar();
            pylab.figure((z + 2)*nc + i + 1)
            pylab.imshow(im[i,0,:,:], vmin = -2, vmax = 2)
            pylab.colorbar();
            pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
