import argparse
import math
import os
import pylab
import sys
import time

BUILD_PATH = os.environ.get('BUILD_PATH') + '/xGadgetron'
SRC_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/pGadgetron'
DATA_PATH = os.environ.get('SRC_PATH') + '/xGadgetron/examples/'
IPT_PATH = os.environ.get('SRC_PATH') + '/ismrmrd-python-tools/ismrmrdtools'

sys.path.append(BUILD_PATH)
sys.path.append(SRC_PATH)
sys.path.append(IPT_PATH)

from pGadgetron import *
from ismrmrdtools import simulation, coils, show

parser = argparse.ArgumentParser(description = \
'''
Test script for CSMs and coil images
''')
parser.add_argument\
('filename', nargs='?', default = 'testdata.h5', \
 help = 'raw data file name (default: testdata.h5)')
args = parser.parse_args()

def csm_sum_range(data):
    u = numpy.conj(data)*data
    u = numpy.sum(u, axis = 0)
    minu = numpy.amin(u)
    maxu = numpy.amax(u)
    return minu, maxu

try:
 
    input_data = MR_Acquisitions(DATA_PATH + args.filename)
    input_data.sort()
    readout = input_data.slice_dimensions()[0]
##    print(input_data.slice_dimensions())

    processed_data = MR_remove_x_oversampling(input_data)
    nx = processed_data.slice_dimensions()[0]
##    print(processed_data.slice_dimensions())

    print(readout, nx)

    print('sorting acquisitions...')
    processed_data.sort()

    cis = MR_CoilImages()

    print('computing coil images...')
##    cis.calculate(input_data)
    cis.calculate(processed_data)

    nz = cis.number()
    print('%d slices' % nz)

    print('computing sensitivity maps...')
    csms = MR_CoilSensitivityMaps()

##    csms.calculate(input_data)
##    csms.calculate(processed_data)

    csms.calculate(cis, method = 'Inati(iter = 1)')
##    csms.calculate(cis, method = '(niter = 10)')

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
        re, im = csms.csm_as_arrays(z)
        zdata = re + 1j*im
##        print(csm_sum_range(zdata))
##        print(minu, maxu)
##    print(minv, maxv)

    nc, m, ny, nx = data.shape
    images = numpy.ndarray((2, ny, nx), dtype = data.dtype)
    allcsms = numpy.ndarray((2*nc, ny, nx), dtype = data.dtype)
    while True:
        s = str(input('enter z-coordinate: '))
        if len(s) < 1:
            break
        z = int(s)
        if z < 0 or z >= nz:
            break
        re, im = cis.coil_image_as_arrays(z)
        coil_data = numpy.squeeze(re + 1j*im)
        (csm, rho) = coils.calculate_csm_inati_iter(coil_data, niter = 10)
##        csm = simulation.generate_birdcage_sensitivities(ny)
##        print(csm_sum_range(csm))

        re, im = csms.csm_as_arrays(z)
        csm_data = numpy.squeeze(re + 1j*im)
        images[0, :, :] = abs(numpy.sum(csm_data * coil_data, axis = 0))
        images[1, :, :] = abs(numpy.sum(csm * coil_data, axis = 0))
        maxv = numpy.amax(images)
        show.imshow(images, tile_shape=(1,2), scale = (0, maxv))

        allcsms[ 0 :   nc, :, :] = abs(csm_data)
        allcsms[nc : 2*nc, :, :] = abs(csm)
        show.imshow(allcsms, tile_shape=(4,4), scale=(0,1))
##        show.imshow(abs(csm), tile_shape=(4,2), scale=(0,1))
##        show.imshow(abs(csm_data), tile_shape=(4,2), scale=(0,1))

##        for i in range(nc):
##            pylab.figure(z*nc + i + 1)
##            pylab.imshow(data[i,0,:,:], vmin = 0, vmax = 1)
####            pylab.figure((z + 1)*nc + i + 1)
####            pylab.imshow(re[i,0,:,:], vmin = -1, vmax = 1)
####            for iy in range(ny):
####                for ix in range(nx):
####                    im[i,0,iy,ix] = math.atan2(im[i,0,iy,ix], re[i,0,iy,ix])
####            pylab.figure((z + 2)*nc + i + 1)
####            pylab.imshow(im[i,0,:,:], vmin = -1, vmax = 1)
##            pylab.show()

except error as err:
    # display error information
    print ('Gadgetron exception occured:\n', err.value)
