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
    csms = MR_CoilSensitivityMaps()
    csm_file = str(input('csm file: '))
    print('reading sensitivity maps...')
    csms.read(csm_file)

    nz = csms.number()
    print('%d slices' % nz)

    maxv = 0
    for z in range(nz):
        data = csms.abs_as_ndarray(z)
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
        data = csms.abs_as_ndarray(z)/maxv
##        data = csms.csm_as_array(z)/maxv
        #shape = data.shape
        zdata = csms.as_ndarray(z)/maxv
        nx, ny, m, nc = csms.map_dimensions()
        re = zdata.real.copy()
        im = zdata.imag.copy()
##        re, im = csms.csm_as_arrays(z)/maxv        
##        shape = re.shape
##        nc = shape[0]
##        ny = shape[1]
##        nx = shape[2]
        print(nx, ny, nc)
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
